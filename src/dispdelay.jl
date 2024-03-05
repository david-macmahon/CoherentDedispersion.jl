const KDM = 4.148808

"""
    dispdelay(flo, fhi, dm)

Return the dispersion delay, in milliseconds, for the given dispersion measure
`dm` from frequency `fhi` to frequency `flo`, which should be given in GHz.
"""
function dispdelay(flo, fhi, dm)
    KDM * dm * (flo^-2 - fhi^-2)
end
