# Bernstein basis

export bernstein_products
"""
Overlap integrals between Bernstein polynomials (which are neither
orthogonal nor normalized)
"""
function bernstein_products(::Type{T}, ::Val{D}, ::Val{P}) where {T,D,P}
    N = D + 1

    scheme = integration_scheme(T, Val(D), Val(P))

    nα = binomial(P + D, D)
    bdb = Array{T}(undef, nα, nα)
    i2α = Array{SVector{N,Int}}(undef, nα)

    αmin = CartesianIndex(ntuple(d -> 0, N))
    αmax = CartesianIndex(ntuple(d -> P, N))

    i = 0
    for αiI in αmin:αmax
        αi = SVector(αiI.I)
        sum(αi) == P || continue
        i += 1

        i2α[i] = αi

        j = 0
        for αjI in αmin:αmax
            αj = SVector(αjI.I)
            sum(αj) == P || continue
            j += 1

            αj >= αi || continue

            f(λ) = bernstein(αi, λ) * bernstein(αj, λ)
            bdb[i, j] = integrate(f, scheme)
        end
        @assert j == nα
    end
    @assert i == nα

    return i2α, Symmetric(bdb)
end

bernstein_products(T, D, P) = bernstein_products(T, Val(D), Val(P))

@generated function integration_scheme(::Type{T}, ::Val{D},
                                       ::Val{order}) where {T,D,order}
    return grundmann_moeller(T, Val(D), order + iseven(order))
end
