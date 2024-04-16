function _copytask(pqin, pqout; id)
    # While non-empty items keep coming in
    while true
        consume!(pqin) do item::NamedTuple
            # If item is empty, return nothing
            isempty(item) && return nothing

            # Send downstream
            produce!(pqout) do data
                copyto!(data, item.data)
                # Synchronize if data needs it (I think copyto! is
                # auto-synchronizing anyway so this is probably not needed)
                coddsynchronize(data)
                return (; data)
            end

            # Return item.data for recycling
            return item.data
        end === nothing && break
    end

    # Send end of input indicator downstream
    produce!(pqout, (;))

    @debug "copytask $id done"
    return nothing
end
