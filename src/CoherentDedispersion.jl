module CoherentDedispersion

using FFTW, LinearAlgebra, Blio, PoolQueues

export KDM, dispdelay

include("dispdelay.jl")
include("coddtask.jl")
include("outputtask.jl")

const CODDInputArrayPQ = PoolQueue{NTuple{2, Array{ComplexF32,2}}, NamedTuple}
const CODDOutputArrayPQ = PoolQueue{NTuple{2, Array{ComplexF32,2}}, NamedTuple}

"""
    coddsynchronize(a)

Potential synchronization point for `a`.  Default method is a no-op, but this
function may be overloaded for specific types of `a` (e.g. a GPU Array type) if
synchronization is required/desired.
"""
function coddsynchronize(a)
end

end #module CoherentDedispersion