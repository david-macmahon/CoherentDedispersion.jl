"""
    dispfreq(dt, fhi, dm)

Return the frequency that is `dt` seconds after `fhi` MHz for the given
dispersion measure `dm`.  Returns `fhi` if `dt` is 0.
"""
function dispfreq(dt, fhi, dm)
    dt == 0 && return fhi
    fhi <= 0 && error("fhi must be >0")
    dm == 0 && error("dm cannot be 0 when dt is non-zero")

    # dt = KDM * dm * (flo^-2 - fhi^-2)
    #
    # flo^-2 = dt / (KDM*dm) + fhi^-2
    #
    #   1        dt        1     dt*fhi^2 + (KDM*dm)
    # ----- = -------- + ----- = -------------------
    # flo^2   (KDM*dm)   fhi^2     KDM * dm * fhi^2
    #
    #           KDM * dm * fhi^2
    # flo^2 = -------------------
    #         dt*fhi^2 + (KDM*dm)
    #
    #            /   KDM * dm * fhi^2  \
    # flo = sqrt(  -------------------  )
    #            \ dt*fhi^2 + (KDM*dm) /
    sqrt(KDM * dm * fhi^2 / (dt*fhi^2 + KDM*dm))
end
