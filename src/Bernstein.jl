"""
Bernstein basis functions for simplices
"""
module Bernstein

using ComputedFieldTypes
using GrundmannMoeller
using LinearAlgebra
using StaticArrays

struct Identity{T}
    data::T
end
Base.getindex(x::Identity) = x.data

include("barycentric.jl")
include("polynomials.jl")
include("basis.jl")

end
