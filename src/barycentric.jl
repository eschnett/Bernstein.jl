# Barycentric coordinates



export cartesian2barycentric
"""
Convert Cartesian to barycentric coordinates

The simplex has one vertex at the origin, the other vertices at unit
distance along the cardianal axes.
"""
function cartesian2barycentric(p::Chain{V, 1, T}) where {V, T}
    D = ndims(V)
    d = 1 - sum(p.v)
    D1 = D+1
    S1 = Signature(D1)
    V1 = SubManifold(S1)
    Chain{V1, 1}(d, p.v...)::Chain{V1, 1, T}
end

export barycentric2cartesian
"""
Convert barycentric to Cartesian coordinates

The simplex has one vertex at the origin, the other vertices at unit
distance along the cardianal axes.
"""
function barycentric2cartesian(λ::Chain{V1, 1, T}) where {V1, T}
    D1 = ndims(V1)
    D = D1-1
    S = Signature(D)
    V = SubManifold(S)
    #TODO <https://github.com/JuliaArrays/StaticArrays.jl/issues/791>
    Chain{V, 1}(SVector{D}([λ.v[i+1] for i in 1:D]))::Chain{V, 1, T}
end



function cartesian2barycentric_setup(s::SVector{N, <:Chain{V, 1, T}}
                                     ) where {N, V, T}
    D = ndims(V)
    @assert N == D+1
    # Algorithm as described on
    # <https://en.wikipedia.org/wiki/Barycentric_coordinate_system>,
    # section "Conversion between barycentric and Cartesian
    # coordinates":
    D1 = D+1
    S1 = Signature(D1)
    V1 = SubManifold(S1)
    extend(x) = Chain{V1, 1}(x.v..., T(1))
    s1 = map(extend, s)
    #TODO <https://github.com/JuliaArrays/StaticArrays.jl/issues/791>
    A = SMatrix{D1, D1}([s1[j][i] for i in 1:D1, j in 1:D1])
    inv(A)::SMatrix{N, N, T}
end

function cartesian2barycentric(invA::SMatrix{N, N, T}, p::Chain{V, 1, T}
                               ) where {N, V, T}
    D = ndims(V)
    @assert N == D+1
    D1 = D+1
    S1 = Signature(D1)
    V1 = SubManifold(S1)
    extend(x) = Chain{V1, 1}(x.v..., T(1))
    b = extend(p).v
    x = invA * b
    Chain{V1, 1}(x)::Chain{V1, 1, T}
end

function cartesian2barycentric(s::SVector{N, <:Chain{V, 1, T}},
                               p::Chain{V, 1, T}
                               ) where {N, V, T}
    invA = cartesian2barycentric_setup(s)
    cartesian2barycentric(invA, p)
end

# # TODO: Use dual numbers for this
# function grad_off_cartesian2barycentric(s::SVector{N, <:Chain{V, 1, T}}
#                                         ) where {N, V, T}
#     invA = cartesian2barycentric_setup(s)
#     D = ndims(V)
#     @assert N == D+1
#     #TODO <https://github.com/JuliaArrays/StaticArrays.jl/issues/791>
#     SMatrix{N, N-1}([invA[i,j] for i in 1:N, j in 1:D]),
#     #TODO <https://github.com/JuliaArrays/StaticArrays.jl/issues/791>
#     SVector{N}([invA[i,D+1] for i in 1:N])
# end

function barycentric2cartesian(s::SVector{N, <:Chain{V, 1, T}},
                               λ::Chain{V1, 1, T}
                               ) where {N, V, V1, T}
    D = ndims(V)
    @assert N == D+1
    D1 = ndims(V1)
    @assert D1 == D+1
    sum(λ[n] * s[n] for n in 1:N)::Chain{V, 1, T}
end
