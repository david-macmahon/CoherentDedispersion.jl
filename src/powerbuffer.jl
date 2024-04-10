# Power buffer
#
# CODD Power buffer:
# CODD integration buffer (mapreducedim! dest): dims = (Nfpc, Nchan, 1, Ntpo)
#     autos::NTuple{2, AbstractArray{Float32,4}}
#     cross::AbstractArray{ComplexF32,4}
# CODD integration buffer output reshaped views: dims = (Nfpc, Nchan, Ntpo)
#     autos::NTuple{2, AbstractArray{Float32,3}}
#     cross::AbstractArray{ComplexF32,3}

#abstract type AbstractCODDPowerBuffer end

struct CODDPowerBuffer{T4<:AbstractArray{Float32,4},
                       Z4<:AbstractArray{ComplexF32,4},
                       T3<:AbstractArray{Float32,3},
                       Z3<:AbstractArray{ComplexF32,3}} #<: AbstractCODDPowerBuffer
    autos4d::NTuple{2,T4}
    cross4d::Z4
    autos::NTuple{2,T3}
    cross::Z3
end

function CODDPowerBuffer(::Type{T}, nfpc, nchan, ntpo) where {T<:AbstractArray}
    autos4d = ntuple(i->T{Float32}(undef, nfpc, nchan, 1, ntpo), 2)
    cross4d = T{ComplexF32}(undef, nfpc, nchan, 1, ntpo)
    autos = reshape.(autos4d, Ref((nfpc, nchan, ntpo)))
    cross = reshape(cross4d, nfpc, nchan, ntpo)
    fill!(CODDPowerBuffer(autos4d, cross4d, autos, cross), 0)
end

function CODDPowerBuffer(nfpc, nchan, ntpo)
    CODDPowerBuffer(Array, nfpc, nchan, ntpo)
end

function Base.fill!(dst::CODDPowerBuffer, x)
    foreach(pol->fill!(pol, x), dst.autos4d)
    fill!(dst.cross4d, x)
    dst
end

function Base.copyto!(dst::CODDPowerBuffer, src::CODDPowerBuffer)
    foreach(((d,s),)->copyto!(d,s), zip(dst.autos4d, src.autos4d))
    foreach(((d,s),)->copyto!(d,s), zip(dst.cross4d, src.cross4d))
    dst
end

function detect!(dst::CODDPowerBuffer, src::CODDVoltageBuffer)
    foreach(zip(dst.autos4d, src.preints)) do (d,s)
        Base.mapreducedim!(abs2, +, d, s)
    end
    bc_conj = Broadcast.broadcasted(conj, src.preints[2])
    bc_cross = Broadcast.broadcasted(*, src.preints[1], bc_conj)
    Base.reducedim!(+, dst.cross4d, Broadcast.instantiate(bc_cross))
    dst
end
