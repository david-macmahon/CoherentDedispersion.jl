const CODDVoltagePQCtor = PoolQueue{CODDVoltageBuffer, NamedTuple}
const CODDPowerPQCtor = PoolQueue{CODDPowerBuffer, NamedTuple}

# PoolQueues are a bit fickle about Channel{T} vs T
function create_plans(cvpq::PoolQueue{<:AbstractChannel{CODDVoltageBuffer}})
    cvb = acquire!(cvpq)
    codd_plan  = plan_fft!(cvb.inputs[1], 1)
    upchan_plan = plan_fft!(cvb.upchans[1], 1);
    recycle!(cvpq, cvb)

    (; codd_plan, upchan_plan)
end

"""
    create_poolqueues_plans(ntpi, nfpc, nint, ntpo, nchan;
                            N=2, use_cuda=CUDA.functional())

Create the PoolQueues and FFT plans for a CoherentDedispersion pipeline.  If
`use_cuda` is `false`, PoolQueues and FFT plans for a CPU-only pipeline will be
created.  If `use_cuda` is true (and CUDA is functional), PoolQueues and FFT
plans for a CPU/GPU pipeline will be created.  `N` is the number of items to
populate each PoolQueue with.

# Extended Help

CPU-only pipeline:

    Input Task (pushes overlapping blocks of RAW data)
    |
    |
     > CPU Voltage PQ
    |
    V
    CODD Task (CODD/upchan/detect)
    |
    |
     > CPU Power PQ
    |
    V
    Output Task (writes Filterbank files)

CPU/GPU Pipeline:

    Input Task (pushes overlapping blocks of RAW data)
    |
    |
     > CPU Voltage PQ
    |
    V
    GPU Copy Task (H2D)
                      |
                      |
                       > GPU Voltage PQ
                      |
                      V
                      CODD Task (CODD/upchan/detect)
                      |
                      |
                       > GPU Power PQ
                      |
                      V
    GPU Copy Task (D2H)
    |
    |
     > CPU Power PQ
    |
    V
    Output Task (writes Filterbank files)
"""
function create_poolqueues_plans(ntpi, nfpc, nint, ntpo, nchan;
                                 N=2, use_cuda=CUDA.functional())
    # cvpq = CPU Voltage PoolQueue
    cvpq = CODDVoltagePQCtor(N) do
        cvb = CODDVoltageBuffer(Array, ntpi, nfpc, nint, ntpo, nchan)
        if use_cuda
            foreach(pol->Mem.pin(pol), cvb.inputs)
        end
        cvb
    end

    # cppq = CPU Power PoolQueue
    cppq = CODDPowerPQCtor(N) do
        cpb = CODDPowerBuffer(Array, nfpc, nchan, ntpo)
        if use_cuda
            foreach(pol->Mem.pin(pol), cpb.autos4d)
            Mem.pin(cpb.cross4d)
        end
        cpb
    end

    if use_cuda
        # gvpq = GPU Voltage PoolQueue
        gvpq = CODDVoltagePQCtor(N) do
            CODDVoltageBuffer(CuArray, ntpi, nfpc, nint, ntpo, nchan)
        end

        # gppq = GPU Power PoolQueue
        gppq = CODDPowerPQCtor(N) do
            CODDPowerBuffer(CuArray, nfpc, nchan, ntpo)
        end

        codd_plan, upchan_plan = create_plans(gvpq)

        return (; cvpq, gvpq, gppq, cppq), codd_plan, upchan_plan
    else
        codd_plan, upchan_plan = create_plans(cvpq)

        return (; cvpq, cppq), codd_plan, upchan_plan
    end
end

function create_tasks(blks, pqs;
                      ntpi, dtpi, fbname, fbheader,
                      f0j, dfj, dm, codd_plan, upchan_plan,
                      use_cuda=CUDA.functional())
    inputtask = errormonitor(
        Threads.@spawn _inputtask(blks, pqs.cvpq;
                                  ntpi, dtpi, fbname, fbheader)
    )

    pqin, pqout = if use_cuda
        pqs.gvpq, pqs.gppq
    else
        pqs.cvpq, pqs.cppq
    end
    coddtask = errormonitor(
        Threads.@spawn _coddtask(pqin, pqout;
                                 f0j, dfj, dm, codd_plan, upchan_plan)
    )

    outputtask = errormonitor(
        Threads.@spawn _outputtask(pqs.cppq)
    )

    tasks = if use_cuda
        htodtask = errormonitor(
            Threads.@spawn _copytask(pqs.cvpq, pqs.gvpq; id="HtoD")
        )

        dtohtask = errormonitor(
            Threads.@spawn _copytask(pqs.gppq, pqs.cppq; id="DtoH")
        )

        (; inputtask, htodtask, coddtask, dtohtask, outputtask)
    else
        (; inputtask, coddtask, outputtask)
    end

    available = Threads.nthreads()
    desired = length(tasks) + 1 # +1 for main thread
    if available < desired
        @warn "only have $available thread(s) for $desired tasks"
    end

    return tasks
end

function create_pipeline(rawfiles, dm;
                         nfpc=1, nint=4, outdir=".", use_cuda=CUDA.functional())
    # Validate use_cuda (in case it is passed explicitly)
    if !CUDA.functional()
        @warn "CUDA is not functional, using CPU-only"
        use_cuda = false
    elseif !use_cuda
        @info "CUDA is functional, but CPU-only requested"
    else
        @info "CUDA is functional and will be used"
    end

    # Load RAW files (reads headers, mmap's data blocks).  Use explicitly typed
    # Vector for datablocks so we don't have to refine later.
    hdrs, blks = GuppiRaw.load(rawfiles; datablocks=Array{Complex{Int8},3}[]);

    # Fixup first header
    CoherentDedispersion.fixup!(hdrs[1])

    npol, = size(blks[1], 1)
    npol == 2 || error("only dual-pol files are supported")

    nchan = size(blks[1], 3)

    # Compute ntime values
    ntpi, ntpo = compute_ntimes(hdrs[1], dm; nfpc, nint)
    # dtpi is number of time samples to advance per block (ntpi-overlap)
    dtpi = nfpc * nint * ntpo

    # Get freq info from header
    obsfreq = hdrs[1][:obsfreq]
    obsbw = hdrs[1][:obsbw]

    # Get parameters for H!
    f0j = obsfreq - obsbw/2
    dfj = obsbw / (hdrs[1].obsnchan / get(hdrs[1], :nants, 1))

    # Create PoolQueues and FFT plans
    pqs,
    codd_plan,
    upchan_plan = create_poolqueues_plans(ntpi, nfpc, nint, ntpo, nchan; use_cuda)

    # Create filterbank filename from name of first rawfile
    rawbase = basename(first(rawfiles))
    fbbase = replace(rawbase, r"\d\d\d\d.raw$"=>"rawcodd.0000.fil")
    fbname = joinpath(outdir, fbbase)

    # Create filterbank header (with values updated based on nfpc and nint)
    # rawspec used machine_id=20 so we might as well use machine_id=21 (it's
    # basically an anachronism at this point).  We also force nifs=4 and then
    # tweak some values to account for the upchannelization and time
    # integration (if any).
    fbheader = Filterbank.Header(hdrs[1]; machine_id=21, nifs=4,
                                          rawdatafile=rawbase)
    fbheader[:nchans] *= nfpc
    fbheader[:foff] /= nfpc
    fbheader[:fch1] -= (nfpc÷2) * fbheader[:foff]
    fbheader[:tsamp] *= nfpc*nint
    fbheader[:nifs] = 4
    # For compatiblity with rawspec, use RA_STR and DEC_STR to get RA/dec.
    # Blio.jl opts not to take on a dependency for the HH:MM:SS.s parsing so it
    # uses the dynamic (i.e. from the telecope's encoders) RA/DEC fields.
    if haskey(hdrs[1], :ra_str)
        fbheader[:src_raj] = hms2ha(hdrs[1][:ra_str])
    end
    if haskey(hdrs[1], :dec_str)
        fbheader[:src_dej] = dms2deg(hdrs[1][:dec_str])
    end

    # TODO create pqs and plans inside create_tasks
    tasks = create_tasks(blks, pqs;
                         ntpi, dtpi, fbname, fbheader,
                         f0j, dfj, dm, codd_plan, upchan_plan, use_cuda)

    (pqs, tasks)
end

function create_pipeline(rawfile::AbstractString, dm;
                         nfpc=1, nint=4, outdir=".", use_cuda=CUDA.functional())
    create_pipeline([rawfile], dm; nfpc, nint, outdir, use_cuda)
end
