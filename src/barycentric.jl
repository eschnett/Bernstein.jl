# Barycentric coordinates

export cartesian2barycentric
"""
    λ = cartesian2barycentric(s, p)

Convert Cartesian to barycentric coordinates.

# Arguments
- `s`: Simplex vertices in Cartesian coordinates. `s` has `N ≤ D + 1`
  vertices in `D` dimensions.
- `p`: Point in Cartesian coordinates

# Result
- `λ`: Point in barycentric coordinates
"""
function cartesian2barycentric(s::SMatrix{N,D,T}, p::SVector{D,T}) where {N,D,T}
    @assert N ≤ D + 1
    # Algorithm as described on
    # <https://en.wikipedia.org/wiki/Barycentric_coordinate_system>,
    # section "Conversion between barycentric and Cartesian
    # coordinates"
    A = SMatrix{N,N}(i == D + 1 ? T(1) : s[j, i] for i in 1:N, j in 1:N)
    b = SVector{D + 1}(p..., T(1))
    return A \ b
end

function cartesian2barycentric(s::SVector{N,SVector{D,T}},
                               p::SVector{D,T}) where {D,N,T}
    return cartesian2barycentric(SMatrix{N,D,T}(s[i][j] for i in 1:N, j in 1:D),
                                 p)
end

export BarycentricSetup
@computed struct BarycentricSetup{N,T}
    invA::fulltype(SMatrix{N,N,T})
end

export cartesian2barycentric_setup
"""
    cartesian2barycentric_setup(s)

Prepare to convert Cartesian to barycentric coordinates.

The returned `setup` structure can be passed to
`cartesian2barycentric` instead of the simplex vertices. This
pre-calculates certain operations and is more efficient.

# Arguments
- `s`: Simplex vertices in Cartesian coordinates. `s` has `N=D+1`
  vertices in `D` dimensions.

# Result
- `setup`
"""
function cartesian2barycentric_setup(s::SMatrix{N,D,T}) where {N,D,T}
    @assert N == D + 1
    A = SMatrix{N,N,T}(j == D + 1 ? T(1) : s[i, j] for j in 1:N, i in 1:N)
    return BarycentricSetup{N,T}(inv(A))
end

function cartesian2barycentric_setup(s::SVector{N,SVector{D,T}}) where {D,N,T}
    return cartesian2barycentric_setup(SMatrix{N,D,T}(s[i][j]
                                                      for i in 1:N, j in 1:D))
end

function cartesian2barycentric(setup::BarycentricSetup{N,T},
                               p::SVector{D,T}) where {N,D,T}
    @assert N == D + 1
    b = SVector{D + 1}(p..., T(1))
    return setup.invA * b
end

export barycentric2cartesian
"""
    p = barycentric2cartesian(s, λ)

Convert barycentric to Cartesian coordinates.

# Arguments
- `s`: Simplex vertices in Cartesian coordinates. `s` has `N ≤ D + 1`
  vertices in `D` dimensions.
- `λ`: Point in barycentric coordinates

# Result
- `p`: Point in Cartesian coordinates
"""
function barycentric2cartesian(s::SMatrix{N,D,T}, λ::SVector{N,T}) where {N,D,T}
    @assert N ≤ D + 1
    return s' * λ
end

function barycentric2cartesian(s::SVector{N,SVector{D,T}},
                               λ::SVector{N,T}) where {N,D,T}
    return barycentric2cartesian(SMatrix{N,D,T}(s[i][j] for i in 1:N, j in 1:D),
                                 λ)
end
