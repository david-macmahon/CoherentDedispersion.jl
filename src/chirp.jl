"""
    function H(f, f0, dm)

Compute a coherent dedispersion complex phase factor at frequency `f+f0`
relative to `f0` for dispersion measure `dm`.  `f` and `f0` should be given in
MHz.
"""
function H(f, f0, dm)
    cispi(2e6 * KDM * dm * f^2 / (f0^2 * (f+f0)))
end

function H(f::Float32, f0::Float32, dm::Float32)
    cispi(2f6 * KDM32 * dm * f^2 / (f0^2 * (f+f0)))
end

"""
    function H(i::Integer, j::Integer, ni, f0j, dfj, dm, dfi=dfj/ni)
    function H(ij::CartesianIndex, ni, f0j, dfj, dm, dfi=dfj/ni)

Compute a coherent dedispersion complex phase factor for a single frequency.
The frequency is determined by a two dimensional indexing scheme that
corresponds to the indexing of a Matrix whose columns contain the Fourier
transform of band limited time domain voltage data.  The columns can be thought
of as "coarse channels" and the elements of each column can be thought of as
"fine channels".  `i` is the zero-based offset of the fine channel.  `j` is the
zero-based offset of the coarse channel.  Alternatively, `ij` is the one-based
`CartesianIndex` corresponding to `(i+1, j+1)`.  `dfi` is the bandwidth of each
fine channel. `ni` is the total number of fine channels.  `f0j` is the reference
frequency for the first coarse channel.  `dfj` is the bandwidth of each coarse
channel.  `dm` is the dispersion measure being dedispersed.

This function assumes that `dfi` is negative and `f0j` refers to the upper edge
of the first coarse channel (or `dfi` is positive and `f0j` refers to the lower
edge of the first coarse channel).  The phase factor is computed relative to the
upper edge of each coarse channel rather than across all coarse channels.  Other
schemes are possible, but they are probably best implemented in a similar but
separate method (rather than as conditional branches within this function).
"""
function H(i::Integer, j::Integer, ni, f0j, dfj, dm, dfi=dfj/ni)
    # Compute FFT shifted value of i.  For top edge f0j and ni == 8, we want ii
    # to go: [4, 5, 6, 7, 0, 1, 2, 3]
    ii = (i + (ni÷2)) % ni

    # Compute f
    f = ii * dfi

    # Compute f0 for coarse channel j
    f0 = f0j + j * dfj

    H(f, f0, dm)
end

function H(ij::CartesianIndex, ni::Integer, f0j::Float32, dfj::Float32, dm::Float32, dfi::Float32=dfj/ni)
    H(ij[1]-1, ij[2]-1, ni, f0j, dfj, dm, dfi)
end

function H(ij::CartesianIndex, ni, f0j, dfj, dm, dfi=dfj/ni)
    H(ij[1]-1, ij[2]-1, Int(ni), Float32(f0j), Float32(dfj), Float32(dm), Float32(dfi))
end

"""
    function H!(m::AbstractMatrix, f0j, dfj, dm; ni=size(m,1), dfi=dfj/ni)

Apply coherent dedispersion complex phase factors to the elements of `m`.  This
essentially multiples each element of `m` by `H(ij, ni, f0j, dfj, dm, dfi)`,
where `ij` is the element's `CartesianIndex`.  `ni` is the size of the first
dimension of `m`.  `dfi` is `dfj/ni`.
"""
function H!(m::AbstractMatrix, f0j, dfj, dm; ni=size(m,1), dfi=dfj/ni)
    m .*= H.(CartesianIndices(m), ni, f0j, dfj, dm, dfi)
end
