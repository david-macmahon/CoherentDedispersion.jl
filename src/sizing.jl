"""
Structure that holds sizing info for a CODDPipeline.  Fields are:

* `ntpi::Int` - Number of time samples per input unit processed
* `nfpc::Int` - Number of fine channels per coarse channel
* `nint::Int` - Number of time samples to integrate
* `ntpo::Int` - Number of time samples output per unit processed
* `nchan::Int` - Number of coarse channels
* `dm::<:AbstractFloat` - Dispersion measure being dedispersed
"""
struct CODDPipelineSize{T<:AbstractFloat}
    ntpi::Int
    nfpc::Int
    nint::Int
    ntpo::Int
    nchan::Int
    dm::T
end

# Construct an invalid CODDPipelineSize instance
CODDPipelineSize() = CODDPipelineSize(0,0,0,0,0,0)

"""
    CODDPipelineSize(grh, dm; [ntpi,] nfpc=1, nint=4) -> CODDPipelineSize

Compute `ntpi` and `ntpo` given dispersion measure `dm` to be coherently
dedispersed and sizing/frequency info in GuppiRaw `grh`.  `ntpi` is the number
of input time samples processed each iteration, `ntpo` is the number of valid
time samples produced each iteration.  `ntpo` will be less than `ntpi` due to
dispersion correction rendering some lower frequency samples at the end of the
input timespan invalid.

`ntpi` defaults to the number of time samples in a single GUPPI RAW block (as
specified by `grh`) but it may be given explicitly as a keyword argument if
desired.  `nfpc` is the number of fine channels per coarse channel to produce.
Pass `nfpc=1` to disable upchannelization.  `nint` is the number of time samples
that will be integrated after detection (i.e. conversion to power).  Due to some
matrix resizing constraints, `nfpc*nint` must divide `ntpi` (and will divide
`ntpo`).  Currently, `ntpi` cannot exceed the number of time samples in a single
GUPPI RAW block, but it may be smaller if desired (though this is not
recommended).
"""
function CODDPipelineSize(grh::GuppiRaw.Header, dm;
                          ntpi=GuppiRaw.getntime(grh), nfpc=1, nint=4)
    # nfpc*nint must divide ntpi (due to reshape tricks in CODDVoltageBuffer)
    ntpi % (nfpc*nint) == 0 || error("nfpc*nint must divide ntpi")

    # Fixup header
    fixup!(grh)

    npol = get(grh, :npol, 1)
    if npol == 1
        @error "only dual-pol files are supported"
        return CODDPipelineSize()
    end

    nants = get(grh, :nants, 1)
    if nants != 1
        @error "only single antennas files are supported"
        return CODDPipelineSize()
    end

    nchan = get(grh, :obsnchan)

    # Compute freq info
    obsbw = grh[:obsbw]
    flo = grh[:obsfreq] - abs(obsbw)/2
    # Some old/broken files don't have chan_bw
    foff = get(grh, :chan_bw, obsbw/nchan)

    # maxdd is dispersion delay across lowest freq coarse channel
    maxdd = dispdelay(flo, flo+abs(foff), dm)

    # overlap is how many time samples correspond to maxdd.  This is the
    # number of samples that must be carried forward to the next
    # workbuf.
    overlap = cld(maxdd, grh[:tbin]) |> Int

    # ntpi is the number of input time samples processed per buffer
    # ntpo is the number of valid output time samples per buffer
    ntpo = ntpi - overlap
    # Express ntpo in integer units of granularity (i.e. nfpc*nint)
    granularity = nfpc * nint
    ntpo = fld(ntpo, granularity)

    # Make sure ntpi is big enough (by making sure ntpo_all is large enough)
    ntpo_all = ntpo * granularity
    if ntpo_all < 0
        @error "dispersion delay exceeds input buffer duration"
        return CODDPipelineSize()
    elseif ntpo_all < ntpi/10
        warn("dispersion delay exceeds 90% input buffer duration")
    end

    return CODDPipelineSize(ntpi, nfpc, nint, ntpo, nchan, dm)
end

"""
    isvalid(cpsz::CODDPipelineSize) -> Bool

Returns true if `cpsz` is valid (i.e. all fields, excluding `dm`, are non-zero).
"""
function isvalid(cpsz::CODDPipelineSize)::Bool
    cpsz.ntpi  != 0 &&
    cpsz.nfpc  != 0 &&
    cpsz.nint  != 0 &&
    cpsz.ntpo  != 0 &&
    cpsz.nchan != 0
end
