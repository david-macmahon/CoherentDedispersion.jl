# Voltage buffer
#
# GUPPI RAW Format Array{Complex{Int8},3} (Npol,Ntpi,Nchan)
#
# CODD Voltage buffer:
# CODD buffer input format NTuple{2,Matrix{Complex{Float32}}} ((Ntpi,Nchan), (Ntpi,Nchan))
# CODD buffer output reshaped view ((Nfpc,Nint,Ntpo,Nchan), (Nfpc,Nint,Ntpo,Nchan))
# CODD buffer pre-integration PermutedDimsArray (1,4,2,3) == (Nfpc, Nchan, Nint, Ntpo)

#import Base: copyto!

#abstract type AbstractCODDVoltageBuffer end

struct CODDVoltageBuffer{T2<:AbstractArray{ComplexF32,2},
                         T4<:AbstractArray{ComplexF32,4}} #<: AbstractCODDVoltageBuffer
    inputs::NTuple{2,T2}
    upchans::NTuple{2,T4}
    preints::NTuple{2,PermutedDimsArray{ComplexF32,4}}
end

function CODDVoltageBuffer(::Type{T}, ntpi, nfpc, nint, ntpo, nchan) where {T<:AbstractArray}
    inputs = ntuple(i->T{ComplexF32}(undef, ntpi, nchan), 2)
    upchans = view.(reshape.(inputs, Ref((nfpc, nint, :, nchan))), :, :, Ref(1:ntpo), :)
    preints = PermutedDimsArray.(upchans, Ref((1,4,2,3)))
    CODDVoltageBuffer(inputs, upchans, preints)
end

function CODDVoltageBuffer(ntpi, nfpc, nint, ntpo, nchan)
    CODDVoltageBuffer(Matrix, ntpi, nfpc, nint, ntpo, nchan)
end

function copyraw!(dst::CODDVoltageBuffer, src::AbstractArray{<:Complex{<:Integer}})
    foreach(((d,s),)->copyto!(d, s), zip(dst.inputs, eachslice(src, dims=1)))
    dst
end

function Base.copyto!(dst::CODDVoltageBuffer, src::CODDVoltageBuffer)
    foreach(((d,s),)->copyto!(d,s), zip(dst.inputs, src.inputs))
    dst
end

function H!(cvb::CODDVoltageBuffer, f0j, dfj, dm; ni=size(cvb.inputs[1],1), dfi=dfj/ni)
    foreach(pol->H!(pol, f0j, dfj, dm; ni, dfi), cvb.inputs)
    cvb
end

function coddsynchronize(cvb::CODDVoltageBuffer)
    coddsynchronize(cvb.inputs[1])
end
