module Speckles

using Reexport
@reexport using JuliaDB
@reexport using UUIDs
@reexport using Logging
@reexport using LaTeXStrings
@reexport using Plots
@reexport using Dates
@reexport using Random
@reexport using StatsBase
@reexport using Distributions
@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using IterTools
@reexport using CSV
@reexport using FFTW


results_directory = "results"

include("Parameters.jl")
include("CountGenerator.jl")
include("Noise.jl")
include("LightSource.jl")
include("Detector.jl")
include("SpecialFunctions.jl")
include("FourierTransform.jl")
include("Strings.jl")
include("FileSave.jl")


#-------------------------------------------------------------------------------

"""
    DenseResults(beamData::IndexedTable, corrData::IndexedTable)

Stores beam intensity and correlation data from Speckle simulation in a dense format.
"""
struct DenseResults
    beamData::IndexedTable
    corrData::IndexedTable
end

#-------------------------------------------------------------------------------
"""
    SparseResults(beamData::NDSparse, corrData::NDSparse)

Stores beam intensity and correlation data from Speckle simulation in a sparse format.
"""
struct SparseResults
    beamData::NDSparse
    corrData::NDSparse
end

struct SpeckleSim
    source::LightSource
    detect::Detector
    bs::Beamsplitter
    duration::Number
    γint::γIntensity
    params::SpeckleParams
end

function SpeckleSim(params::SpeckleParams)
    source = LightSource(
                         params.n,
                         params.Em,
                         params.νm,
                         params.σ,
                         params.fγ
                        )
    detect = Detector(
                      params.deadtime,
                      params.resolution,
                      params.jitter,
                      params.efficiency,
                      params.darkcounts
                     )
    bs = Beamsplitter(1,1)
    nb = nbar(params.duration,source)*bs.t
    γint = γIntensity(nb,params.duration,detect.resolution,source)

    return SpeckleSim(source,  detect, bs, params.duration, γint, params)
end

export SpeckleSim
#-------------------------------------------------------------------------------
function run(sim::SpeckleSim)
    # id = uuid4()
    # # caulcate the max index of the time integration window
    # nwindow = length(sim.γint)÷convert(Int, ceil(sim.params.duration/sim.params.window))
    
    # # generate correlation times
    # τ = sim.params.resolution*collect(0:(nwindow-1))

    # # create a 5 e-fold low cut in τ to get rid of the non-periodic part of g2
    # τCutLow = sqrt(5)/sim.params.σ

    # # index of τCutLow
    # nτCutLow = convert(Int, ceil(τCutLow/sim.params.resolution))

    # # get the frequencies of a Fourier transform over τ with the calculated cuts
    # freqs = fftFreq(sim.params.resolution,τ,(nτCutLow,lastindex(τ)))

    return "Hello, World!"



end
export run

#-------------------------------------------------------------------------------
end
