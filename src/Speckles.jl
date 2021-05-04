module Speckles

using Reexport
@reexport using Random
@reexport using StatsBase
@reexport using Distributions
@reexport using DataFrames
@reexport using LinearAlgebra
@reexport using IterTools
@reexport using CSV
@reexport using FFTW

include("CountGenerator.jl")
include("SpeckleNoise.jl")
include("ClassicalEM.jl")
include("SpeckleFunctions.jl")
include("SpeckleFT.jl")

# struct SpeckleParams
    # params::eFieldParams
    # seed::Integer
    # tres::Number
    # tmax::Number
    # windowCounts::Integer
    # totalCounts::Integer
# end

# export SpeckleParams

struct Beamsplitter
    r::Number
    t::Number
    function Beamsplitter(r::Number,t::Number)
        beamNorm = sqrt(r^2+t^2)
        new(r/beamNorm,t/beamNorm)
    end
end

export Beamsplitter

end
