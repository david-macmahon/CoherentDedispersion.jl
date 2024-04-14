module CoherentDedispersion

using FFTW, LinearAlgebra, Blio, PoolQueues, CUDA, CUDA.CUFFT
using RadioInterferometry # For guppifixup.jl
using ProgressBars

export KDM, KDM32, dispdelay, dispfreq
export CODDVoltageBuffer, copyraw!
export CODDPowerBuffer, detect!
export H!
export compute_ntimes
export create_pipeline

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
include("copytask.jl")
include("outputtask.jl")

include("guppifixup.jl")
include("coddpipeline.jl")

"""
    coddsynchronize(a)

Potential synchronization point for `a`.  Default method is a no-op, but this
function may be overloaded for specific types of `a` (e.g. a GPU Array type) if
synchronization is required/desired.
"""
function coddsynchronize(_)
end

function coddsynchronize(_::CuArray)
    synchronize()
end

end #module CoherentDedispersion
