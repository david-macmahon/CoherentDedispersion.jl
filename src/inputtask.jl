function _inputtask(data, pqout; ntpi, ntpo, fbname, fbheader)
    # Get total number of time samples from data
    ntime = size(data, 2)
    
    # Loop thrrough input time samples, sending ntpi time samples each
    # iteration, but only advancing ntpo time samples each iteration.
    for t in 1:ntpo:(ntime-ntpi+1)
        # Send data downstream
        produce!(pqout) do cvb
            copyraw!(cvb, @view data[:,t:t+ntpi-1])
            coddsynchronize(cvb) # Synchronize if cvb needs it
            return (; fbname, fbheader, cvb)
        end
    end

    # Send end of input indicator downstream
    produce!(pqout, (;))

    @info "inputtask done"
    return nothing
end
