# Voltage buffer
#
# GUPPI RAW Format Array{Complex{Int8},3} (Npol,Ntpi,Nchan)
#
# CODD Voltage buffer:
# CODD buffer input format NTuple{2,Matrix{Complex{Float32}}} ((Ntpi,Nchan), (Ntpi,Nchan))
# CODD buffer output reshaped view ((Nfpc,Nint,:,Nchan), (Nfpc,Nint,Ntpi÷(Nfpc*Nint),Nchan))
# CODD buffer pre-integration SubArray of PermutedDimsArray (1,4,2,3) == (Nfpc, Nchan, Nint, Ntpo)

struct CODDVoltageBuffer{T2<:AbstractArray{ComplexF32,2},
                         T4<:AbstractArray{ComplexF32,4}}
    inputs::NTuple{2,T2}
    upchans::NTuple{2,T4}
    preints::NTuple{2,SubArray{ComplexF32,4}}
end

function CODDVoltageBuffer(::Type{T}, ntpi, nfpc, nint, ntpo, nchan) where {T<:AbstractArray}
    inputs = ntuple(i->T{ComplexF32}(undef, ntpi, nchan), 2)
    # upchans are reshaped views of inputs that we will upchannelize along the
    # first dimention (of size nfpc).  We will end up upchannelizing the entire
    # input buffer rather than just the the first nfpc*nint*ntpo time samples,
    # but this allows CUDA.CUFFT to work with the data in situ rather than
    # having to copy to another buffer so it's a worthwhile trade-off (IMHO, I
    # have not benchmarked it).
    upchans = reshape.(inputs, Ref((nfpc, nint, :, nchan)))
    preints = view.(PermutedDimsArray.(upchans, Ref((1,4,2,3))), :, :, :, Ref(1:ntpo))
    CODDVoltageBuffer(inputs, upchans, preints)
end

function CODDVoltageBuffer(ntpi, nfpc, nint, ntpo, nchan)
    CODDVoltageBuffer(Matrix, ntpi, nfpc, nint, ntpo, nchan)
end

function copyraw!(dst::CODDVoltageBuffer,
                  blks::AbstractVector{<:AbstractArray{<:Complex{<:Integer}}},
                  tidx, ntime=size(dst.inputs[1], 1))
    # 11111111 22222222 33333333    blk_idx
    # 12345678 12345678 12345678    blk_tidx
    #               ^               tidx=14
    b0, n = divrem(tidx-1, ntime) # b0=1, n=5 for tidx=14 and ntime=8

    nblks = length(blks)
    if b0+1 > nblks
        for d in dst.inputs
            fill!(d, zero(eltype(d)))
        end
        return dst
    end

    # Copy from first block
    src = blks[b0+1]
    axchan = axes(dst.inputs[1],2)
    srcidxs = CartesianIndices((n+1:ntime, axchan))
    dstidxs = CartesianIndices((1:ntime-n, axchan))
    for (d,s) in zip(dst.inputs, eachslice(src, dims=1))
        copyto!(d, dstidxs, s, srcidxs)
    end

    # Copy from second block, if needed
    if n > 0
        dstidxs = CartesianIndices((ntime-n+1:ntime, axchan))
        if b0+2 <= nblks
            src = blks[b0+2]
            srcidxs = CartesianIndices((1:n, axchan))
            for (d,s) in zip(dst.inputs, eachslice(src, dims=1))
                copyto!(d, dstidxs, s, srcidxs)
            end
        else
            # Emulate dummy block of zeros at the end of blks
            for d in dst.inputs
                d[dstidxs] .= 0
            end
        end
    end

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
