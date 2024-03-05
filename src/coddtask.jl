"""
    chirp(k::Integer, l::Integer, dm, fhi, foff, Tms, N) -> ComplexF32
    chirp(kl::CartesianIndex, dm, fhi, foff, Tms, N)

Generate chirp phase factor used in `chirp!`.  `k` and `l` are zero-based
offsets.  `kl` is a one-based `CartesianIndex`.  `dm` is the dispersion measure,
`fhi` is the highest frequency in the highest frequency channel.  `foff` is the
frequency offset between coarse channels (and the bandwidth of a single coarse
channel).  `Tms` is the duration in milliseconds of the time span being chirped.
`N` is the number of fine channels per coarse channel.
"""
function chirp(k::Integer, l::Integer, dm, fhi, foff, Tms, N)
    # Compute fhil (fhi for coarse channel l)
    fhil = fhi + l * foff
    # Compute dfk (fine channel offset from fhil)
    dfk = (foff / N) * ((k + N÷2) % N)
    # Compute fkl, the frequency for (k,l)
    fkl = fhil + dfk
    # Compute dispersion delay from fhil to fkl
    dd = dispdelay(fkl, fhil, dm)
    # Return chirp phase factor
    cispi(-2 * fkl * dd / Tms)
end

function chirp(kl::CartesianIndex, dm, fhi, foff, Tms, N)
    chirp(kl[1]-1, kl[2]-1, dm, fhi, foff, Tms, N)
end

"""
    chirp!(workbuf::AbstractMatrix, dm, fhi, foff, Tms)

Apply a dedispersion "chirp" to `workbuf` based on the dispersion measure `dm`,
the highest frequency (of the highest frequency channel) `fhi`, the frequency
offset between coarse channels `foff`, and the total de-chirp time span `Tms`.
`fhi` and `foff` must be given in GHz.  `Tms` must be given in milliseconds.
"""
function chirp!(workbuf::AbstractMatrix, dm, fhi, foff, Tms)
    N = size(workbuf, 1)
    workbuf .*= chirp.(CartesianIndices(workbuf), dm, fhi, foff, Tms, N)
end

function _coddtask(pqin, pqout; workbufs, overlap_granularity=4)
    old_blocks = nothing

    # Create FFT plans
    fplan = plan_fft!(workbufs[1], 1)
    bplan = plan_ifft!(workbufs[1], 1)

    ntime, nchan = size(workbufs[1])
    next = ntime+1
    local dm
    local fhi
    local foff
    local Tms

    # While non-empty items keep coming in
    while true
        consume!(pqin) do item::NamedTuple
            # If item is empty, return nothing
            isempty(item) && return nothing

            # Get dm, header, and block from item
            dm = item.dm
            rawname = item.rawname
            header = item.header
            new_blocks = item.data
            @info "coddtask got dm $dm pktidx $(header[:pktidx]) data $(size(new_blocks[1])) $(new_blocks[1][1])"

            # Recompute freq info etc
            fhi = (header[:obsfreq] + abs(header[:obsbw])/2) / 1e3
            foff = header[:chan_bw] / 1e3
            Tms = header[:tbin] * ntime * 1e3
            maxdd = dispdelay(fhi-nchan*abs(foff), fhi-(nchan-1)*abs(foff), dm)
            overlap = cld(maxdd/1e3, header[:tbin]) |> Int
            # Round overlap up to next multiple of overlap_granularity
            overlap = overlap_granularity * cld(overlap, overlap_granularity)

            # Process old block
            @info "> next $next"
            while next <= ntime
                # Copy from old_blocks
                ntimeold = ntime - next + 1
                rdest = CartesianIndices((1:ntimeold, 1:nchan))
                rsrc  = CartesianIndices((next:ntime, 1:nchan))
                @assert size(rdest) == size(rsrc) "overlap size mismatch"
                copyto!.(workbufs, Ref(rdest), old_blocks, Ref(rsrc))

                # Copy from new_blocks
                rdest = CartesianIndices((ntimeold+1:ntime, 1:nchan))
                rsrc = CartesianIndices((1:next-1, 1:nchan))
                @assert size(rdest) == size(rsrc) "non-overlap size mismatch"
                copyto!.(workbufs, Ref(rdest), new_blocks, Ref(rsrc))

                # Forward FFT block
                mul!.(workbufs, Ref(fplan), workbufs)
                # Apply chirp
                chirp!.(workbufs, dm, fhi, foff, Tms)
                # Backward FFT block
                mul!.(workbufs, Ref(bplan), workbufs)

                # Output workbuf
                produce!(pqout) do bufouts
                    copyto!.(bufouts, workbufs)
                    coddsynchronize.(bufouts) # Synchronize if bufouts need it
                    return (; dm, overlap, rawname, header, data=bufouts)
                end

                # Advance next
                next += ntime - overlap
                @info "- next $next"
            end

            # Out with the old, in with the new
            next -= ntime
            @info "< next $next"
            recyclable = old_blocks
            old_blocks = new_blocks

            # Return `recyclable`.  The first time through, `recyclable` (was
            # `old_blocks`) is nothing, but we can't recycle nothing in that
            # case because that would cause us to break out of the while loop.
            # PoolQueues should support returning Bool from this function where
            # `true` would signify "nothing to recycle, but don't stop and
            # `false` would signify "nothing to recycle and stop" (like
            # `nothing` signifies now).  But until that happens, we acquire! a
            # dummy block from `pqin` and recycle! it.
            if recyclable === nothing
                recyclable = acquire!(pqin)
            end

            return recyclable
        end === nothing && break
    end

    # TODO copy/process remaining portion of old_blocks

    # Recycle old block, produce empty NamedTuple to signify end of input
    recycle!(pqin, old_blocks)
    produce!(pqout, (;))

    @info "coddtask done"
    return nothing
end
    