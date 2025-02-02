using Pkg
Pkg.activate(".")

using LinearAlgebra
using Plots


mutable struct Body
	position::Vector[Float64]
	velocity::Vector[Float64]
	mass::Float64
	force::Vector[Float64]
	radius::Float64
end


const G = 1.0
const θ = 0.5
const dt = 0.1


N = 1000

i = range(1, N)
