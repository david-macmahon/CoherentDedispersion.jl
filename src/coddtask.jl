function _coddtask(pqin, pqout;
    workbufs::NTuple{2,<:AbstractArray{<:Complex,2}},
    overlap_granularity=4
)
    old_blocks = nothing

    # Create FFT plans
    fplan = plan_fft!(workbufs[1], 1)
    bplan = plan_ifft!(workbufs[1], 1)

    ntime, nchan = size(workbufs[1])
    next = ntime+1
    local dm
    local fhi
    local foff

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
            fhi = (header[:obsfreq] + abs(header[:obsbw])/2)
            flo = (header[:obsfreq] - abs(header[:obsbw])/2)
            foff = header[:chan_bw]
            # maxdd is dispersion delay across lowest freq coarse channel
            maxdd = dispdelay(flo+abs(foff), flo, dm)
            # overlap is how many time samples correspond to maxdd.  This is the
            # number of samples that must be carried forward to the next
            # workbuf.
            overlap = cld(maxdd, header[:tbin]) |> Int
            # Round overlap up to next multiple of overlap_granularity.
            # TODO Actually, we should probably round (ntime-overlap) down to
            # the previous multiple of overlap_granularity.
            overlap = overlap_granularity * cld(overlap, overlap_granularity)

            # Parameters for `H!`
            ni = ntime
            dfj = foff
            dfi = dfj/ni
            f0j = fhi

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
                # Apply phase factors
                H!.(workbufs, dm, fhi, foff, Tms)
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
    