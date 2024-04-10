"""
    compute_ntimes(dm, header; [ntpi,] nfpc=1, nint=4) -> (ntpi, ntpo)

Compute `ntpi` and `ntpo` given dispersion measure `dm` to be coherently
dedispersed and sizing/frequency info in GuppiRaw `header`.  `ntpi` is the
number of input time samples processed each iteration, `ntpo` is the number of
valid time samples produced each iteration.  `ntpo` will be less than `ntpi` due
to dispersion correction rendering some lower frequency samples at the end of
the input timespan invalid.

`ntpi` defaults to the number of time samples in a single GUPPI RAW block (as
specified by `header`) but it may be given explicitly as a keyword argument if
desired.  `nfpc` is the number of fine channels per coarse channel to produce.
Pass `nfpc=1` to disable upchannelization.  `nint` is the number of time samples
that will be integrated after detection (i.e. conversino to power).  Due to some
matrix resizing constraints, `nfpc*nint` must divide `ntpi` (and will divide
`ntpo`).
"""
function compute_ntimes(header, dm;
                        ntpi=GuppiRaw.getntime(header), nfpc=1, nint=4)
    # nfpc*nint must divide ntpi (due to reshape tricks in CODDVoltageBuffer)
    ntpi % (nfpc*nint) == 0 || error("nfpc*nint must divide ntpi")

    # Compute freq info
    flo = (header[:obsfreq] - abs(header[:obsbw])/2)
    # Some old/broken files don't have chan_bw
    foff = get(header, :chan_bw, header[:obsbw]/header[:obsnchan])

    # maxdd is dispersion delay across lowest freq coarse channel
    maxdd = dispdelay(flo, flo+abs(foff), dm)

    # overlap is how many time samples correspond to maxdd.  This is the
    # number of samples that must be carried forward to the next
    # workbuf.
    overlap = cld(maxdd, header[:tbin]) |> Int

    # ntpi is the number of input time samples processed per buffer
    # ntpo is the number of valid output time samples per buffer
    ntpo = ntpi - overlap
    # Express ntpo in integer units of granularity (i.e. nfpc*nint)
    granularity = nfpc * nint
    ntpo = fld(ntpo, granularity)

    # Make sure ntpi is big enough (by making sure ntpo_all is large enough)
    ntpo_all = ntpo * granularity
    if ntpo_all < 0
        error("dispersion delay exceeds input buffer duration")
    elseif ntpo_all < ntpi/10
        warn("dispersion delay exceeds 90% input buffer duration")
    end

    return ntpi, ntpo
end
