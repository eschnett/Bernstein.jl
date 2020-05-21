# Bernstein polynomials

export bernstein
"""
Evaluate Bernstein polynomial

    s: Vertices of simplex
    α: multi-index describing polynomial
    x: location where Bernstein polynomial is evaluated

The order of approximation is given implicity by `sum(α)`.
"""
bernstein

function bernstein(s::SVector{N, <: Chain{V, 1, T}},
                   α::SVector{N, UInt},
                   λ::Chain{V1, 1, T}
                   ) where {N, V, V1, T}
    D = ndims(V)
    @assert N == D+1
    D1 = ndims(V1)
    @assert D1 == D+1
    n = sum(α)

    # TODO: Handle factorials better
    T(factorial(n)) * prod(λ[i]^α[i] / T(factorial(α[i])) for i in 1:N)::T
end

function bernstein(s::SVector{N, <: Chain{V, 1, T}},
                   α::SVector{N, UInt},
                   x::Chain{V, 1, T}
                   ) where {N, V, T}
    D = ndims(V)
    @assert N == D+1

    λ = cartesian2barycentric(s, x)
    bernstein(s, α, λ)
end

function bernstein(s::SVector{N, <: Chain{V, 1, T}},
                   α::SVector{N, Int},
                   x::Chain{V1, 1, T}
                   ) where {N, V, V1, T}
    bernstein(s, UInt.(α), x)
end



# function bernstein_barycentric_coefficients(s::SVector{N, <: Chain{V, 1, T}},
#                                             α::SVector{N, UInt}
#                                             ) where {N, V, T}
#     D = ndims(V)
#     @assert N == D+1
#     n = sum(α)
# 
#     # TODO: Handle factorials better
#     nfact = T(factorial(n))
#     #TODO <https://github.com/JuliaArrays/StaticArrays.jl/issues/791>
#     SVector{N, Tuple{T, T}}([(nfact / factorial(α[i]), α[i]) for i in 1:N])
# end
