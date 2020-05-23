"""
Bernstein basis functions for simplices
"""
module Bernstein

using ComputedFieldTypes
using LinearAlgebra
using SimplexQuad
using StaticArrays



include("barycentric.jl")
include("polynomials.jl")
include("basis.jl")

end
