function _inputtask(data, pqout; ntpi, dtpi, fbname, fbheader)
    # Sanity check dtpi vs ntpi
    dtpi > ntpi && error("dtpi must not exceed ntpi")

    # Get total number of time samples from data
    ntime = size(data, 2)
    
    startidxs = 1:dtpi:(ntime-ntpi+1)
    # Loop thrrough input time samples, sending ntpi time samples each
    # iteration, but only advancing dtpi time samples each iteration.
    for (i,t) in enumerate(startidxs)
        @info "processing overlapped block $i starting at time index $t"
        # Send data downstream
        produce!(pqout) do cvb
            copyraw!(cvb, @view data[:,t:t+ntpi-1,:])
            coddsynchronize(cvb) # Synchronize if cvb needs it
            return (; fbname, fbheader, cvb)
        end
    end

    # Send end of input indicator downstream
    produce!(pqout, (;))

    @info "inputtask done"
    return nothing
end
