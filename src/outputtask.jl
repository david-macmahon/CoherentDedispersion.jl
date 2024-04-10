function _outputtask(pqin)
    local fbio
    local auto_fftshift
    local cross_fftshift

    """
        c2r(c) -> r

    Creates a ReinterpretArray for c with a leading dimensino of 2 for re/im.
    If c has size (S1, S2, ...) then r will have size (2, S1, S2, ...).
    """
    function c2r(c::AbstractArray{<:Complex})
        reinterpret(reshape, real(eltype(c)), c)
    end

    # While non-empty items keep coming in
    while true
        consume!(pqin) do item::NamedTuple
            # If item is empty, return nothing
            isempty(item) && return nothing

            # Get CODDPowerBuffer from item
            cpb = item.cpb

            # Open output file and allocate fftshift buffers, if needed
            if !@isdefined fbio
                fbio = open(item.fbname, "w")
                write(fbio, item.fbheader)

                auto_fftshift = similar(cpb.autos[1])
                cross_fftshift = similar(cpb.cross)
                reim_fftshift = reinterpret(reshape,
                    real(eltype(cross_fftshift)), cross_fftshift
                )
            end

            for pol in cbp.autos
                if size(pol, 1) > 1
                    fftshift!(auto_fftshift, pol, 1)
                    write(fbio, auto_fftshift)
                else
                    write(fbio, pol)
                end
            end

            if size(cpb.cross, 1) > 1
                # Upchannelization, do fftshift and write
                fftshift!(cross_fftshift, cpb.cross, 1)
                for reim in eachslice(c2r(cross_fftshift), dims=1)
                    write(fbio, reim)
                end
            else
                # No upchannelization, write directly from cpb.cross
                for reim in eachslice(c2r(cpb.cross), dims=1)
                    write(fbio, reim)
                end
            end

            # Return cpb for recycling
            return cpb
        end === nothing && break
    end

    # Send end of input indicator downstream
    produce!(pqout, (;))

    @info "coddtask done"
    return nothing
end
