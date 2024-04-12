const CODDVoltagePQ = PoolQueue{CODDVoltageBuffer, NamedTuple}
const CODDPowerPQ = PoolQueue{CODDPowerBuffer, NamedTuple}

"""
    create_poolqueues(::Type{Array}, ntpi, nfpc, nint, ntpo, nchan)

Create the PoolQueues for a CoherentDedispersion pipeline.

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
"""
function create_poolqueues(T::Type{<:AbstractArray},
                           ntpi, nfpc, nint, ntpo, nchan; N=2)
    cvpq = CODDVoltagePQ(N) do
        CODDVoltageBuffer(T, ntpi, nfpc, nint, ntpo, nchan)
    end

    cppq = CODDPowerPQ(N) do 
        CODDPowerBuffer(T, nfpc, nchan, ntpo)
    end

    (; cvpq, cppq)
end

function create_poolqueues(ntpi, nfpc, nint, ntpo, nchan)
    create_poolqueues(Array, ntpi, nfpc, nint, ntpo, nchan)
end

# PoolQueues are a bit fickle about Channel{T} vs T
function create_plans(cvpq::PoolQueue{<:AbstractChannel{CODDVoltageBuffer}})
    cvb = acquire!(cvpq)
    codd_plan  = plan_fft!(cvb.inputs[1], 1)
    upchan_plan = plan_fft!(cvb.upchans[1], 1);
    recycle!(cvpq, cvb)

    (; codd_plan, upchan_plan)
end


function create_tasks(data, pqs;
                      ntpi, dtpi, fbname, fbheader,
                      f0j, dfj, dm, codd_plan, upchan_plan)
    inputtask = errormonitor(
        Threads.@spawn _inputtask(data, pqs.cvpq;
                                  ntpi, dtpi, fbname, fbheader)
    )
    coddtask = errormonitor(
        Threads.@spawn _coddtask(pqs.cvpq, pqs.cppq;
                                 f0j, dfj, dm, codd_plan, upchan_plan)
    )
    outputtask = errormonitor(
        Threads.@spawn _outputtask(pqs.cppq)
    )

    (; inputtask, coddtask, outputtask)
end

function create_pipeline(T::Type{<:AbstractArray}, rawfiles, dm;
                         nfpc=1, nint=4, outdir=".")
    # Load data files (reads headers, mmap's data blocks).  Use explicitly typed
    # Vector for datablocks so we don't have to refine later.
    hdrs, blks = GuppiRaw.load(rawfiles; datablocks=Array{Complex{Int8},3}[]);

    # Fixup first header
    CoherentDedispersion.fixup!(hdrs[1])

    npol, = size(blks[1], 1)
    npol == 2 || error("only dual-pol files are supported")

    nchan = size(blks[1], 3)

    # Concatenate all blocks into a BlockArray "view" along the time (i.e.
    # second) dimension
    data = mortar(reshape(blks, (1,:,1)));

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

    # Create PoolQueues
    pqs = create_poolqueues(T, ntpi, nfpc, nint, ntpo, nchan)

    # Create FFT plans
    @show typeof(pqs.cvpq)
    codd_plan, upchan_plan = create_plans(pqs.cvpq)

    # Create filterbank filename from name of first rawfile
    rawname = first(rawfiles)
    fbbase = replace(basename(rawname), r"\d\d\d\d.raw$"=>"rawcodd.0000.fil")
    fbname = joinpath(outdir, fbbase)

    # Create filterbank header (with values updated based on nfpc and nint)
    fbheader = Filterbank.Header(hdrs[1])
    fbheader[:nchans] *= nfpc
    fbheader[:foff] /= nfpc
    fbheader[:fch1] -= (nfpc÷2) * fbheader[:foff]
    fbheader[:tsamp] *= nfpc*nint
    fbheader[:tstart] += fbheader[:tsamp] / (24*60*60) / 2
    fbheader[:nifs] = 4

    tasks = create_tasks(data, pqs;
                         ntpi, dtpi, fbname, fbheader,
                         f0j, dfj, dm, codd_plan, upchan_plan)

    (pqs, tasks)
end

function create_pipeline(T::Type{<:AbstractArray}, rawfile::AbstractString, dm; nfpc=1, nint=4)
    create_pipeline(T, [rawfile], dm; nfpc, nint)
end

function create_pipeline(rawfiles, dm; nfpc=1, nint=4)
    create_pipeline(Array, rawfiles, dm; nfpc, nint)
end

function create_pipeline(rawfile::AbstractString, dm; nfpc=1, nint=4)
    create_pipeline(Array, [rawfile], dm; nfpc, nint)
end
