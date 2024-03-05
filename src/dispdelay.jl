
"""
    dispdelay(flo, fhi, dm)

Return the dispersion delay, in seconds, for the given dispersion measure `dm`
from frequency `fhi` to frequency `flo`, both of which should be given in MHz.
"""
function dispdelay(flo, fhi, dm)
    KDM * dm * (flo^-2 - fhi^-2)
end
