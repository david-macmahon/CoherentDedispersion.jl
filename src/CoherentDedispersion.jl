module CoherentDedispersion

using FFTW, LinearAlgebra, Blio, PoolQueues

export KDM, KDM32, dispdelay, dispfreq
export CODDVoltageBuffer, copyraw!
export CODDPowerBuffer, detect!
export H!
export compute_ntimes

"""
`KDM` is the *dispersion constant* in units of `MHz² pc⁻¹ cm³ s` expressed as a
`Float64`.
"""
const KDM = 4.148808e3

"""
`KDM32` is the *dispersion constant* in units of `MHz² pc⁻¹ cm³ s` expressed as
a `Float32`.
"""
const KDM32 = 4.148808f3

include("chirp.jl")
include("dispdelay.jl")
include("dispfreq.jl")
include("voltagebuffer.jl")
include("powerbuffer.jl")
include("sizing.jl")

include("inputtask.jl")
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