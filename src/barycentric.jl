# Barycentric coordinates



export cartesian2barycentric
"""
Convert Cartesian to barycentric coordinates

The simplex has one vertex at the origin, the other vertices at unit
distance along the cardianal axes.
"""
function cartesian2barycentric(p::SVector{D, T}) where {D, T}
    N = D+1
    d = 1 - sum(p)
    SVector{N}(d, p...)
end

export barycentric2cartesian
"""
Convert barycentric to Cartesian coordinates

The simplex has one vertex at the origin, the other vertices at unit
distance along the cardianal axes.
"""
function barycentric2cartesian(λ::SVector{N, T}) where {N, T}
    D = N-1
    SVector{D, T}(λ[a+1] for a in 1:D)
end



function cartesian2barycentric(s::SMatrix{D, N, T}, p::SVector{D, T}
                               ) where {D, N, T}
    @assert N == D+1
    # Algorithm as described on
    # <https://en.wikipedia.org/wiki/Barycentric_coordinate_system>,
    # section "Conversion between barycentric and Cartesian
    # coordinates":
    A = SMatrix{N, N}(i == D+1 ? T(1) : s[i,j] for i in 1:N, j in 1:N)
    b = SVector{D+1}(p..., T(1))
    A \ b
end

function cartesian2barycentric(s::SVector{N, SVector{D, T}}, p::SVector{D, T}
                               ) where {D, N, T}
    cartesian2barycentric(SMatrix{D, N}(s[j][i] for i in 1:D, j in 1:N), p)
end

export BarycentricSetup
@computed struct BarycentricSetup{N, T}
    invA::fulltype(SMatrix{N, N, T})
end

export cartesian2barycentric_setup
"""
Prepare to convert Cartesian to barycentric coordinates

The returned `setup` structure can be passed to
`cartesian2barycentric` instead of the vertices. This pre-calculates
certain operations and is more efficient.
"""
function cartesian2barycentric_setup(s::SMatrix{D, N, T}) where {D, N, T}
    @assert N == D+1
    A = SMatrix{N, N, T}(i == D+1 ? T(1) : s[i,j] for i in 1:N, j in 1:N)
    BarycentricSetup{N, T}(inv(A))
end

function cartesian2barycentric_setup(s::SVector{N, SVector{D, T}}
                               ) where {D, N, T}
    cartesian2barycentric_setup(SMatrix{D, N, T}(
        s[j][i] for i in 1:D, j in 1:N))
end

function cartesian2barycentric(setup::BarycentricSetup{N, T}, p::SVector{D, T}
                               ) where {N, D, T}
    @assert N == D+1
    b = SVector{D+1}(p..., T(1))
    setup.invA * b
end



# # TODO: Use dual numbers for this
# function grad_off_cartesian2barycentric(s::SVector{N, <:Chain{V, 1, T}}
#                                         ) where {N, V, T}
#     invA = cartesian2barycentric_setup(s)
#     D = ndims(V)
#     @assert N == D+1
#     SMatrix{N, N-1}(invA[i,j] for i in 1:N, j in 1:D],
#     SVector{N}(invA[i,D+1] for i in 1:N)
# end

function barycentric2cartesian(s::SMatrix{D, N, T}, λ::SVector{N, T}
                               ) where {D, N, T}
    @assert N == D+1
    s * λ
end
