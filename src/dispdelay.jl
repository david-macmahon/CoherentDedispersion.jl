const KDM = 4.148808

"""
    dispdelay(dm, flo, fhi)

Return the dispersion delay, in milliseconds, for the given dispersion measure
`dm` from frequency `fhi` to frequency `flo`, which should be given in GHz.
"""
function dispdelay(dm, flo, fhi)
    KDM * dm * (flo^-2 - fhi^-2)
end
