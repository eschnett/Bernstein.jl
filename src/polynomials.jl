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

function bernstein(α::SVector{N,UInt}, λ::SVector{N,T}) where {N,T}
    n = sum(α)

    # TODO: Handle factorials better
    T(factorial(n)) * prod(λ[i]^α[i] / T(factorial(α[i])) for i = 1:N)::T
end

bernstein(α::SVector{N,<:Integer}, λ::SVector{N}) where {N} = bernstein(UInt.(α), λ)



function bernstein(s::SMatrix{D,N,T}, α::SVector{N,UInt}, x::SVector{D,T}) where {D,N,T}
    @assert N == D + 1

    λ = cartesian2barycentric(s, x)
    bernstein(α, λ)
end

function bernstein(
    s::SMatrix{D,N,T},
    α::SVector{N,<:Integer},
    x::SVector{D,T},
) where {D,N,T}
    bernstein(s, UInt.(α), x)
end



export bernstein_setup
bernstein_setup(s) = cartesian2barycentric_setup(s)

function bernstein(
    setup::BarycentricSetup{N,T},
    α::SVector{N,UInt},
    x::SVector{D,T},
) where {N,D,T}
    @assert N == D + 1

    λ = cartesian2barycentric(setup, x)
    bernstein(α, λ)
end

function bernstein(
    setup::BarycentricSetup{N,T},
    α::SVector{N,<:Integer},
    x::SVector{D,T},
) where {D,N,T}
    bernstein(setup, UInt.(α), x)
end



# function bernstein_barycentric_coefficients(s::SMatrix{D, N, T}
#                                             α::SVector{N, UInt}
#                                             ) where {D, N, T}
#     @assert N == D+1
#     n = sum(α)
# 
#     # TODO: Handle factorials better
#     nfact = T(factorial(n))
#     SVector{N, Tuple{T, T}}((nfact / factorial(α[i]), α[i]) for i in 1:N)
# end
