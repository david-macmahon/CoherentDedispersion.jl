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

function CODDPowerBuffer(::Type{T}, sz::CODDPipelineSize) where {T<:AbstractArray}
    CODDPowerBuffer(T, sz.nfpc, sz.nchan, sz.ntpo)
end

function CODDPowerBuffer(nfpc, nchan, ntpo)
    CODDPowerBuffer(Array, nfpc, nchan, ntpo)
end

function CODDPowerBuffer(sz::CODDPipelineSize)
    CODDPowerBuffer(Array, sz)
end

function Base.fill!(dst::CODDPowerBuffer, x)
    foreach(pol->fill!(pol, x), dst.autos4d)
    fill!(dst.cross4d, x)
    dst
end

function Base.copyto!(dst::CODDPowerBuffer, src::CODDPowerBuffer)
    foreach(((d,s),)->copyto!(d,s), zip(dst.autos4d, src.autos4d))
    copyto!(dst.cross4d, src.cross4d)
    dst
end

"""
    sumdiff!(a, b) -> (a, b)

For equally sized Arrays `a` and `b`, modify `a` in place to be `a+b` and modify
`b` in place to be `a-b`.
"""
function sumdiff!(a::AbstractArray, b::AbstractArray)
    for i in eachindex(a,b)
        @inbounds a[i], b[i] = a[i]+b[i], a[i]-b[i]
    end
end

function sumdiff_kernel(a, b)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x

    if i <= length(a) == length(b)
        a[i], b[i] = (a[i] + b[i]), (a[i] - b[i])
    end

    return nothing
end

function sumdiff!(a::CuArray, b::CuArray;
                  threads=min(1024, length(a)),
                  blocks=cld(length(a), threads))
    @cuda blocks threads sumdiff_kernel(a, b)
    a, b
end

function detect!(dst::CODDPowerBuffer, src::CODDVoltageBuffer;
                 dostokes::Bool=true, doconj::Bool=false, doscale::Bool=true)
    # Compute autos (with optional scaling by 127^2 to match rawspec)
    scaling = doscale ? 127f0^2 : 1f0
    foreach(zip(dst.autos4d, src.preints)) do (d,s)
        Base.mapreducedim!(x->abs2(x)/scaling, +, d, s)
    end

    # If stokes is requested, compute I and Q
    dostokes && sumdiff!(dst.autos4d[1], dst.autos4d[2])

    # If inputs are NOT conjugated (conj==false), compute `conj.(pol2).*pol1`.
    # If inputs ARE conjugated (conj==true), compute `conj.(pol2).*pol1`.
    bc_conj = Broadcast.broadcasted(conj, src.preints[2-doconj])
    bc_cross = Broadcast.broadcasted(*, src.preints[1+doconj], bc_conj)
    bc_scale = Broadcast.broadcasted(/, bc_cross, scaling)
    Base.reducedim!(+, dst.cross4d, Broadcast.instantiate(bc_scale))
    dst
end

function coddsynchronize(cpb::CODDPowerBuffer)
    coddsynchronize(cpb.autos4d[1])
end
