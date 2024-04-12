function _coddtask(pqin, pqout; f0j, dfj, dm, codd_plan, upchan_plan)
    # While non-empty items keep coming in
    while true
        consume!(pqin) do item::NamedTuple
            # If item is empty, return nothing
            isempty(item) && return nothing

            # Get fields from item
            fbname = item.fbname
            fbheader = item.fbheader
            cvb = item.data

            # Coherently dedisperse
            foreach(pol->mul!(pol, codd_plan, pol), cvb.inputs)
            H!(cvb, f0j, dfj, dm)
            foreach(pol->mul!(pol, inv(codd_plan), pol), cvb.inputs)

            # Upchannelize (if nchans > 1)
            if size(cvb.preints[1], 1) > 1
                foreach(pol->mul!(pol, upchan_plan, pol), cvb.upchans)
            end

            # Send downstream
            produce!(pqout) do cpb
                fill!(cpb, 0)
                detect!(cpb, cvb)
                coddsynchronize(cpb) # Synchronize if cpb needs it
                return (; fbname, fbheader, data=cpb)
            end

            # Return cvb for recycling
            return cvb
        end === nothing && break
    end

    # Send end of input indicator downstream
    produce!(pqout, (;))

    @info "coddtask done"
    return nothing
end
