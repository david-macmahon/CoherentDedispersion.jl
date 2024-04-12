"""
    c2ri(c) -> ri

Creates a ReinterpretArray for `c` with a leading dimension of 2 for re/im.
If `c` has size (S1, S2, ...) then `ri` will have size (2, S1, S2, ...) with
`real.(c) == ri[1,:,:]` and `imag.(c) == ri[2,:,:]`.
"""
function c2ri(c::AbstractArray{<:Complex})
    reinterpret(reshape, real(eltype(c)), c)
end

function _outputtask(pqin)
    local fbname
    local fbio
    local writebuf
    local a11_writebufvw
    local a22_writebufvw
    local r12_writebufvw
    local i12_writebufvw

    # While non-empty items keep coming in
    while true
        consume!(pqin) do item::NamedTuple
            # If item is empty, return nothing
            isempty(item) && return nothing

            # Get CODDPowerBuffer from item
            cpb = item.cpb
            nfpc, nchan, ntpo = size(cpb.autos[1])

            # Open output file and allocate write buffer, if needed
            if !@isdefined fbname
                fbname = item.fbname
                fbio = open(fbname, "w")
                write(fbio, item.fbheader)

                writebuf = similar(cpb.autos[1], nfpc, nchan, 4*ntpo)
                a11_writebufvw = @view writebuf[:,:,1:4:end]
                a22_writebufvw = @view writebuf[:,:,2:4:end]
                r12_writebufvw = @view writebuf[:,:,3:4:end]
                i12_writebufvw = @view writebuf[:,:,4:4:end]
            end

            # Interleave spectra from autos and cross into writebuf
            fftshift!(a11_writebufvw, cpb.autos[1], 1)
            fftshift!(a22_writebufvw, cpb.autos[2], 1)
            cross_reim = c2ri(cpb.cross)
            fftshift!(r12_writebufvw, @view(cross_reim[1,:,:,:]), 1)
            fftshift!(i12_writebufvw, @view(cross_reim[2,:,:,:]), 1)

            # Write writebuf to Filterbank file
            write(fbio, writebuf)

            # Return cpb for recycling
            return cpb
        end === nothing && break
    end

    @info "coddtask done"

    return fbname
end
