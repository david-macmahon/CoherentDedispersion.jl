function _inputtask(blks, pqout; ntpi, dtpi, progress=false)
    # Sanity check dtpi vs ntpi
    dtpi > ntpi && error("dtpi must not exceed ntpi")

    # Get total number of time samples over all blocks
    ntime = sum(b->size(b,2), blks)
    
    startidxs = 1:dtpi:(ntime-ntpi+1)
    # Loop through input time samples, sending ntpi time samples each
    # iteration, but only advancing dtpi time samples each iteration.
    pbar = progress ? ProgressBar : identity
    for t in pbar(startidxs)
        # Send data downstream
        produce!(pqout) do cvb
            copyraw!(cvb, blks, t)
            return (; data=cvb)
        end
    end

    # Send end of input indicator downstream
    produce!(pqout, (;))

    @debug "inputtask done"
    return nothing
end
