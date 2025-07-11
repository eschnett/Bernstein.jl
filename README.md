# Bernstein Polynomial Basis

* [GitHub](https://github.com/eschnett/Bernstein.jl): Source code repository
* [![GitHub CI](https://github.com/eschnett/Bernstein.jl/workflows/CI/badge.svg)](https://github.com/eschnett/Bernstein.jl/actions)

The [Bernstein
polynomials](https://en.wikipedia.org/wiki/Bernstein_polynomial) form
a basis for polynomials living on simplices. This package calculates
the Bernstein polynomials for simplices of arbitrary dimension.

This package also provides conversion functions between Cartesian and
[barycentric
coordinates](https://en.wikipedia.org/wiki/Barycentric_coordinate_system).

## Examples

### Convert between Cartesian and barycentric coordinates

```julia
julia> using StaticArrays

julia> using Bernstein

julia> # Define a triangle
       s = rand(SMatrix{3,2})
3×2 SArray{Tuple{3,2},Float64,2,6} with indices SOneTo(3)×SOneTo(2):
 0.814346  0.297149
 0.519781  0.620776
 0.345743  0.733385

julia> # Choose a point
       p = rand(SVector{2})
2-element SArray{Tuple{2},Float64,1,2} with indices SOneTo(2):
 0.1483582649665245
 0.3923504628863179

julia> # Convert to barycentric coordinates
       λ = cartesian2barycentric(s, p)
3-element SArray{Tuple{3},Float64,1,3} with indices SOneTo(3):
   3.523588891718682
 -10.621534570787583
   8.097945679068905

julia> # Convert back
       q = barycentric2cartesian(s, λ)
2-element SArray{Tuple{2},Float64,1,2} with indices SOneTo(2):
 0.14835826496652604
 0.3923504628863199

julia> p ≈ q
true
```

You can also pass the simplex vertices as a vector of vectors
`SVector{N, SVector{D, T}}` instead of a matrix.

If you convert many Cartesian to barycentric coordinates, then part of
the transformation can be pre-calculated to increase efficiency. Call
`cartesian2barycentric_setup` for this:

```julia
julia> using StaticArrays

julia> using Bernstein

julia> # Define a triangle
       s = rand(SMatrix{3,2})
3×2 SArray{Tuple{3,2},Float64,2,6} with indices SOneTo(3)×SOneTo(2):
 0.814346  0.297149
 0.519781  0.620776
 0.345743  0.733385

julia> # Choose a point
       p = rand(SVector{2})
2-element SArray{Tuple{2},Float64,1,2} with indices SOneTo(2):
 0.1483582649665245
 0.3923504628863179

julia> # Pre-calculate part of the transformation
       setup = cartesian2barycentric_setup(s);

julia> # Convert to barycentric coordinates
       λ = cartesian2barycentric(setup, p)
3-element SArray{Tuple{3},Float64,1,3} with indices SOneTo(3):
   3.523588891718682
 -10.621534570787583
   8.097945679068905
```

### Evaluate Bernstein polynomials

You can evaluate Bernstein polynomials from barycentric coordinates or
from Cartesian coordinates:

```julia
julia> using StaticArrays

julia> using Bernstein

julia> # Define a triangle
       s = rand(SMatrix{3,2})
3×2 SArray{Tuple{3,2},Float64,2,6} with indices SOneTo(3)×SOneTo(2):
 0.814346  0.297149
 0.519781  0.620776
 0.345743  0.733385

julia> # Choose a point
       p = rand(SVector{2})
2-element SArray{Tuple{2},Float64,1,2} with indices SOneTo(2):
 0.1483582649665245
 0.3923504628863179

julia> # Convert to barycentric coordinates
       λ = cartesian2barycentric(s, p)
3-element SArray{Tuple{3},Float64,1,3} with indices SOneTo(3):
   3.523588891718682
 -10.621534570787583
   8.097945679068905

julia> # Choose polynomial index and order
       # (The order is the sum of all coefficients)
       α = SVector(2,0,0);

julia> bernstein(α, λ)
3.523588891718682

julia> bernstein(s, α, p)
3.523588891718682

julia> bernstein(setup, α, p)
3.523588891718682
```

### 1D Bernstein polynomials

For 1D polynomials are defined in the interval [0,1].
Polynomial of order `n` and index `ν`,

```math
B_{ν,n}(x) = \binom{n}{ν} x^ν (1-x)^{n-ν}
```
in [Wikipedia conventions](https://en.wikipedia.org/wiki/Bernstein_polynomial), can be called as:

```julia
using Bernstein, Bernstein.StaticArrays

function bernstein_1d(ν::Int, n::Int, x::Number)
    λ = SVector(1.0 - x, x)
    α = SVector(n - ν, ν)
    return bernstein(α, λ)
end
```

## References

- Douglas N. Arnold, Richard S. Falk, Ragnar Winther, "Geometric
  decompositions and local bases for spaces of finite element
  differential forms", 10.1016/j.cma.2008.12.017, [arXiv:0806.1255
  [math.NA]](https://arxiv.org/abs/0806.1255).

- Tom Lyche, Karl Scherer, "On the p-norm condition number of the
  multivariate triangular Bernstein basis", [DOI:
  10.1016/S0377-0427(00)00383-6](https://doi.org/10.1016/S0377-0427(00)00383-6).
