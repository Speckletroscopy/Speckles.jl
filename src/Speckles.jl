module Speckles

using Reexport
@reexport using JuliaDB
@reexport using Logging
@reexport using LaTeXStrings
@reexport using Plots
@reexport using Dates
@reexport using Random
@reexport using StatsBase
@reexport using Distributions
@reexport using DataFrames
@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using IterTools
@reexport using CSV
@reexport using FFTW

include("CountGenerator.jl")
include("Noise.jl")
include("LightSource.jl")
include("Detector.jl")
include("SpecialFunctions.jl")
include("FourierTransform.jl")
include("Strings.jl")
include("FileSave.jl")

struct SpeckleSim
    source::LightSource
    detect::Detector
    bs::Beamsplitter
    duration::Number
    γint::γIntensity
end

#-------------------------------------------------------------------------------
"""
    DenseResults(beamData::IndexedTable, corrData::IndexedTable)

Stores beam intensity and correlation data from Speckle simulation in a dense format.
"""
struct DenseResults
    beamData::IndexedTable
    corrData::IndexedTable
end

"""
    SparseResults(beamData::NDSparse, corrData::NDSparse)

Stores beam intensity and correlation data from Speckle simulation in a sparse format.
"""
struct SparseResults
    beamData::NDSparse
    corrData::NDSparse
end

"""

    struct SpeckleParams

Stores all relevant parameters to create a simulation

--- LightSource variables ---
n::UInt = Number of atoms
νm::Vector{Float64}    = Line frequencies in GHz
Em::Vector{Float64}    = Relative line magnitudes
σ::Float64             = Doppler broadening in GHz 
fγ::Float64            = Mean photon cout rate in GHz
--- Detector variables ---
deadtime::Float64      = Detector deadtime in nanoseconds
resolution::Float64    = Detector resolutinon in nanoseconds
jitter::Float64        = Detector jitter in nanoseconds
efficiency::Float64    = Detector efficiency
darkcounts::Float64    = Detector dark count rate in GHz
--- Simulation parameter variables ---
duration::Float64      = Duration of readout in nanoseconds
window::Float64        = Duration of time averaging in nanoseconds
repeat::UInt           = Number of times to repeat correlation measurement
reinstance::Bool       = Randomize frequencies and phases after each repeat
"""
struct SpeckleParams{T<:Float64}
    # LightSource variables
    n::UInt # Number of atoms
    νm::Vector{T} # Line frequencies in GHz
    Em::Vector{T} # Relative line magnitudes
    σ::T # Doppler broadening in GHz 
    fγ::T # Mean photon cout rate in GHz
    # Detector variables
    deadtime::T # Detector deadtime in nanoseconds
    resolution::T # Detector resolutinon in nanoseconds
    jitter::T # Detector jitter in nanoseconds
    efficiency::T # Detector efficiency
    darkcounts::T # Detector dark count rate in GHz
    # Simulation parameter variables
    duration::T # Duration of readout in nanoseconds
    window::T # Duration of time averaging in nanoseconds
    repeat::UInt # T of times to repeat correlation measurement
    reinstance::Bool # Randomize frequencies and phases after each repeat
    function SpeckleParams{T}(n::UInt, νm::Vector{T}, Em::Vector{T}, σ::T, 
                              fγ::T, deadtime::T, resolution::T, jitter::T,
                              efficiency::T, darkcounts::T, duration::T,
                              window::T, repeat::UInt, reinstance::Bool) where T<:Float64
        @assert length(νm) == length(Em) "Line magnitudes and frequencies must match lengths"
        new(n, νm, Em, σ, fγ, deadtime, resolution, jitter, efficiency,
            darkcounts, duration, window, repeat, reinstance)
    end
end

"""
    function SpeckleParams(n::Number, νm::Vector{Number}, Em::Vector{Number}, σ::Number, 
                              fγ::Number, deadtime::Number, resolution::Number, jitter::Number,
                              efficiency::Number, darkcounts::Number, duration::Number,
                              window::Number, repeat::Number, reinstance::Bool)

Outer contructor for the SpeckleParams datatype that allows for more general type input.
"""
function SpeckleParams(n::Number, νm::Vector{Number}, Em::Vector{Number}, σ::Number, 
                          fγ::Number, deadtime::Number, resolution::Number, jitter::Number,
                          efficiency::Number, darkcounts::Number, duration::Number,
                          window::Number, repeat::Number, reinstance::Bool)

    SpeckleParams( convert(UInt, n),
                  convert(Vector{Float64},νm),
                  convert(Vector{Float64},Em),
                  convert(Float64,σ),
                  convert(Float64,fγ),
                  convert(Float64,deadtime),
                  convert(Float64,resolution),
                  convert(Float64,jitter),
                  convert(Float64,efficiency),
                  convert(Float64,darkcounts),
                  convert(Float64,duration),
                  convert(Float64,window),
                  convert(UInt,repeat),
                  reinstance
                 )
end

"""
    function SpeckleParamsVector(params::Dict)

Takes a dictionary of proper keywords with values or vectors of values for a
SpeckleParam struct and splits it into all possible valid combinations of
those values. Returns a vector of SpeckleParams.
"""
function SpeckleParamsVector(params::Dict)
    paramVec = paramVector(params)
    out = map(p->SpeckleParams(p[:n],
                               p[:νm],
                               p[:Em],
                               p[:σ],
                               p[:fγ],
                               p[:deadtime],
                               p[:resolution],
                               p[:jitter],
                               p[:efficiency],
                               p[:darkcounts],
                               p[:duration],
                               p[:window],
                               p[:repeat],
                               p[:reinstance]
                              ),
                              paramVec
                             )

    return out
end

#-------------------------------------------------------------------------------
function keyReplace!(params::Dict)
    # ones
    if params[:Em] == "ones"
        params[:Em] = ones(length(params[:νm]))
    end

    # Counts from poisson error
    if params[:fγ] == "shot1%"
        params[:fγ] = 2*1e4/params[:resolution] # multiply by 2 so error in each beam is ~1%
    elseif params[:fγ] == "shot10%"
        params[:fγ] = 2*1e2/params[:resolution] # multiply by 2 so error in each beam is ~10%
    elseif params[:fγ] == "shot50%"
        params[:fγ] = 2*4/params[:resolution] # multiply by 2 so error in each beam is ~50%
    end

    # correlation window
    if params[:timeint] == "halfwindow"
       params[:timeint] = params[:duration]/2 
    end
    return params
end

function keyReplace(params::Dict)
    newdict = deepcopy(params)
    # ones
    if newdict[:Em] == "ones"
        newdict[:Em] = ones(length(newdict[:νm]))
    end

    # Counts from poisson error
    if newdict[:fγ] == "shot1%"
        newdict[:fγ] = 2*1e4/newdict[:resolution] # multiply by 2 so error in each beam is ~1%
    elseif newdict[:fγ] == "shot10%"
        newdict[:fγ] = 2*1e2/newdict[:resolution] # multiply by 2 so error in each beam is ~10%
    elseif newdict[:fγ] == "shot50%"
        newdict[:fγ] = 2*4/newdict[:resolution] # multiply by 2 so error in each beam is ~50%
    end

    # correlation window
    if newdict[:timeint] == "halfwindow"
       newdict[:timeint] = newdict[:duration]/2 
    end
    return newdict
end
#-------------------------------------------------------------------------------

"""
    paramVector(params::Dict)

Splits multiple parameter dictionary into a vector of individual parameter dictionaries.
"""
function paramVector(params::Dict)
    k = collect(keys(params))
    v = collect(Iterators.product(collect(values(params))...))
    return map(vals->Dict(collect(zip(k,vals))),ivec(v))
end

export paramVector

#-------------------------------------------------------------------------------
end
