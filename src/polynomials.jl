# Bernstein polynomials

export bernstein
"""
    bernstein(α, λ)
    bernstein(s, α, x)

Evaluate Bernstein polynomial

The order of approximation is given implicity by `sum(α)`.

# Arguments
- `s`: Vertices of simplex
- `α`: multi-index describing polynomial
- `x`: location where Bernstein polynomial is evaluated in Cartesian
  coordinates
- `λ`: location where Bernstein polynomial is evaluated in barycentric
  coordinates
"""
bernstein

function bernstein(α::SVector{N,UInt}, λ::SVector{N,T}) where {N,T}
    n = sum(α)

    # TODO: Handle factorials better
    return T(factorial(n)) *
           prod(λ[i]^α[i] / T(factorial(α[i])) for i in 1:N)::T
end

function bernstein(α::SVector{N,<:Integer}, λ::SVector{N}) where {N}
    return bernstein(UInt.(α), λ)
end

function bernstein(s::SMatrix{N,D,T}, α::SVector{N,UInt},
                   x::SVector{D,T}) where {N,D,T}
    @assert N == D + 1

    λ = cartesian2barycentric(s, x)
    return bernstein(α, λ)
end

function bernstein(s::SMatrix{N,D,T}, α::SVector{N,<:Integer},
                   x::SVector{D,T}) where {D,N,T}
    return bernstein(s, UInt.(α), x)
end

export bernstein_setup
bernstein_setup(s) = cartesian2barycentric_setup(s)

function bernstein(setup::BarycentricSetup{N,T}, α::SVector{N,UInt},
                   x::SVector{D,T}) where {N,D,T}
    @assert N == D + 1

    λ = cartesian2barycentric(setup, x)
    return bernstein(α, λ)
end

function bernstein(setup::BarycentricSetup{N,T}, α::SVector{N,<:Integer},
                   x::SVector{D,T}) where {D,N,T}
    return bernstein(setup, UInt.(α), x)
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
