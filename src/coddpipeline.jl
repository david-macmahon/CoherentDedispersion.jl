const CODDVoltagePQCtor = PoolQueue{CODDVoltageBuffer, NamedTuple}
const CODDPowerPQCtor = PoolQueue{CODDPowerBuffer, NamedTuple}

struct CODDPipelineCPU{CV<:AbstractChannel{<:CODDVoltageBuffer},
                       CP<:AbstractChannel{<:CODDPowerBuffer},
                       P2<:AbstractFFTs.Plan, P4<:AbstractFFTs.Plan}
    cvpq::PoolQueue{CV, Channel{NamedTuple}}
    cppq::PoolQueue{CP, Channel{NamedTuple}}
    codd_plan::P2
    upchan_plan::P4
    cpsz::CODDPipelineSize
end

struct CODDPipelineGPU{CV<:AbstractChannel{<:CODDVoltageBuffer},
                       GV<:AbstractChannel{<:CODDVoltageBuffer},
                       CP<:AbstractChannel{<:CODDPowerBuffer},
                       GP<:AbstractChannel{<:CODDPowerBuffer},
                       P2<:AbstractFFTs.Plan, P4<:AbstractFFTs.Plan}
    cvpq::PoolQueue{CV, Channel{NamedTuple}}
    gvpq::PoolQueue{GV, Channel{NamedTuple}}
    gppq::PoolQueue{GP, Channel{NamedTuple}}
    cppq::PoolQueue{CP, Channel{NamedTuple}}
    codd_plan::P2
    upchan_plan::P4
    cpsz::CODDPipelineSize
end

# PoolQueues are a bit fickle about Channel{T} vs T
function create_plans(cvpq::PoolQueue{<:AbstractChannel{CODDVoltageBuffer}})
    cvb = acquire!(cvpq)
    codd_plan  = plan_fft!(cvb.inputs[1], 1)
    upchan_plan = plan_fft!(cvb.upchans[1], 1);
    recycle!(cvpq, cvb)

    (; codd_plan, upchan_plan)
end

# CPU-only pipeline
function _create_pipeline(cpsz::CODDPipelineSize, _use_cuda::Val{false}; N=2)
    # cvpq = CPU Voltage PoolQueue
    cvpq = CODDVoltagePQCtor(N) do
        CODDVoltageBuffer(Array, cpsz)
    end

    codd_plan, upchan_plan = create_plans(cvpq)

    # cppq = CPU Power PoolQueue
    cppq = CODDPowerPQCtor(N) do
        CODDPowerBuffer(Array, cpsz)
    end

    CODDPipelineCPU(cvpq, cppq, codd_plan, upchan_plan, cpsz)
end

# CPU/GPU pipeline
function _create_pipeline(cpsz::CODDPipelineSize, _use_cuda::Val{true}; N=2)
    # cvpq = CPU Voltage PoolQueue
    cvpq = CODDVoltagePQCtor(N) do
        cvb = CODDVoltageBuffer(Array, cpsz)
        for pol in cvb.inputs
            Mem.pin(pol)
        end
        cvb
    end

    # gvpq = GPU Voltage PoolQueue
    gvpq = CODDVoltagePQCtor(N) do
        CODDVoltageBuffer(CuArray, cpsz)
    end

    codd_plan, upchan_plan = create_plans(gvpq)

    # gppq = GPU Power PoolQueue
    gppq = CODDPowerPQCtor(N) do
        CODDPowerBuffer(CuArray, cpsz)
    end

    # cppq = CPU Power PoolQueue
    cppq = CODDPowerPQCtor(N) do
        cpb = CODDPowerBuffer(Array, cpsz)
        for pol in cpb.autos4d
            Mem.pin(pol)
        end
        Mem.pin(cpb.cross4d)
        cpb
    end

    CODDPipelineGPU(cvpq, gvpq, gppq, cppq, codd_plan, upchan_plan, cpsz)
end

"""
    create_pipeline(grh::GuppiRaw.Header, dm;
                    nfpc=1, nint=4, N=2, use_cuda=CUDA.functional())

Create the PoolQueues and FFT plans for a CoherentDedispersion pipeline based on
sizing info `grh`, `dm` (the dispersion measure), `nfpc` (the number of fine
channels per coarse channel), and `nint` (the number of time samples to
integrate).  If `use_cuda` is `false`, PoolQueues and FFT plans for a
CPU-only pipeline will be created.  If `use_cuda` is true (and CUDA is
functional), PoolQueues and FFT plans for a CPU/GPU pipeline will be created.
`N` is the number of items to populate each PoolQueue with.

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
function create_pipeline(grh::GuppiRaw.Header, dm;
                         nfpc=1, nint=4, N=2, use_cuda=CUDA.functional())
    # Validate use_cuda (in case it is passed explicitly)
    if !CUDA.functional()
        @warn "CUDA is not functional, using CPU-only"
        use_cuda = false
    elseif !use_cuda
        @info "CUDA is functional, but CPU-only requested"
    else
        @info "CUDA is functional and will be used"
    end

    cpsz = CODDPipelineSize(grh, dm; nfpc, nint)

    if !isvalid(cpsz)
        error("invalid sizing info detected")
    end

    # Create actual pipeline (i.e. PoolQueues and FFT plans)
    _create_pipeline(cpsz, Val(use_cuda); N)
end

function create_pipeline(rawfile::AbstractString, dm;
                         nfpc=1, nint=4, N=2, use_cuda=CUDA.functional())
    grh = read(rawfile, GuppiRaw.Header)
    create_pipeline(grh, dm; nfpc, nint, N, use_cuda)
end

function create_pipeline(rawfiles::AbstractVector, dm;
                         nfpc=1, nint=4, N=2, use_cuda=CUDA.functional())
    create_pipeline(first(rawfiles), dm; nfpc, nint, N, use_cuda)
end

function _create_tasks(pipeline::CODDPipelineCPU, blks;
                      fbname, fbheader, f0j, dfj, progress=false)
    dm   = pipeline.cpsz.dm
    ntpi = pipeline.cpsz.ntpi
    nfpc = pipeline.cpsz.nfpc
    nint = pipeline.cpsz.nint
    ntpo = pipeline.cpsz.ntpo
    dtpi = nfpc * nint * ntpo

    inputtask = errormonitor(
        Threads.@spawn _inputtask(blks, pipeline.cvpq; ntpi, dtpi, progress)
    )

    coddtask = errormonitor(
        Threads.@spawn _coddtask(pipeline.cvpq, pipeline.cppq; f0j, dfj, dm,
                                 pipeline.codd_plan, pipeline.upchan_plan)
    )

    outputtask = errormonitor(
        Threads.@spawn _outputtask(pipeline.cppq; fbname, fbheader)
    )

    (; inputtask, coddtask, outputtask)
end

function _create_tasks(pipeline::CODDPipelineGPU, blks;
                      fbname, fbheader, f0j, dfj, progress=false)
    dm   = pipeline.cpsz.dm
    ntpi = pipeline.cpsz.ntpi
    nfpc = pipeline.cpsz.nfpc
    nint = pipeline.cpsz.nint
    ntpo = pipeline.cpsz.ntpo
    dtpi = nfpc * nint * ntpo

    inputtask = errormonitor(
        Threads.@spawn _inputtask(blks, pipeline.cvpq; ntpi, dtpi, progress)
    )

    htodtask = errormonitor(
        Threads.@spawn _copytask(pipeline.cvpq, pipeline.gvpq; id="HtoD")
    )

    coddtask = errormonitor(
        Threads.@spawn _coddtask(pipeline.gvpq, pipeline.gppq; f0j, dfj, dm,
                                 pipeline.codd_plan, pipeline.upchan_plan)
    )

    dtohtask = errormonitor(
        Threads.@spawn _copytask(pipeline.gppq, pipeline.cppq; id="DtoH")
    )

    outputtask = errormonitor(
        Threads.@spawn _outputtask(pipeline.cppq; fbname, fbheader)
    )

    (; inputtask, htodtask, coddtask, dtohtask, outputtask)
end

function start_pipeline(pipeline, blks::AbstractVector{<:AbstractArray};
                        fbname, fbheader, f0j, dfj, progress=false)
    tasks = _create_tasks(pipeline, blks; fbname, fbheader, f0j, dfj, progress)

    available = Threads.nthreads()
    desired = length(tasks) + 1 # +1 for main thread
    if available < desired
        @warn "only have $available thread(s) for $desired tasks" maxlog=1
    end

    tasks
end

function start_pipeline(pipeline, rawfiles::AbstractVector{<:AbstractString};
                        outdir=".", progress=false)
    dm   = pipeline.cpsz.dm
    nfpc = pipeline.cpsz.nfpc
    nint = pipeline.cpsz.nint

    # Load RAW files (reads headers, mmap's data blocks), but only if sizing
    # matches.  For now we require all fields to be equal, but this may be
    # loosened in the future to allow `dm` to vary (downward).
    hdrs, blks = GuppiRaw.load(rawfiles) do grh
        pipeline.cpsz === CODDPipelineSize(grh, dm; nfpc, nint)
    end

    # If hdrs is empty, return empty task list indicating problem (e.g. size
    # mismatch)
    isempty(hdrs) && return (;)

    # Fixup first header
    fixup!(hdrs[1])

    # Get freq info from header
    obsfreq = hdrs[1][:obsfreq]
    obsbw = hdrs[1][:obsbw]

    # Get parameters for H!
    f0j = obsfreq - obsbw/2
    dfj = obsbw / (hdrs[1].obsnchan / get(hdrs[1], :nants, 1))

    # TODO Make function to get (fbname, fbheader)
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

    # TODO Return blks as well so that caller can finalize them upon completion?
    start_pipeline(pipeline, blks; fbname, fbheader, f0j, dfj, progress)
end

# TODO Return empty String rather than nothing when tasks in empty?
function run_pipeline(pipeline, blks::AbstractVector{<:AbstractArray};
                      fbname, fbheader, f0j, dfj, progress=true)
    tasks = start_pipeline(pipeline, blks; fbname, fbheader, f0j, dfj, progress)
    # TODO Finalize blks after completion?
    isempty(tasks) ? nothing : fetch(last(tasks))
end

# TODO Return empty String rather than nothing when tasks in empty?
function run_pipeline(pipeline, rawfiles::AbstractVector{<:AbstractString};
                      outdir=".", progress=true)
    tasks = start_pipeline(pipeline, rawfiles; outdir, progress)
    # TODO Finalize blks after completion?
    isempty(tasks) ? nothing : fetch(last(tasks))
end
