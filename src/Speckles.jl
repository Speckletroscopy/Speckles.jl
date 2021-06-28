module Speckles

using Reexport
@reexport using JuliaDB
@reexport using UUIDs
@reexport using Dates
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

# location where all results are stored
results_directory = "results"

include("CountGenerator.jl")
include("Noise.jl")
include("LightSource.jl")
include("Detector.jl")
include("FourierTransform.jl")
include("Parameters.jl")
include("SimulationFlow.jl")
include("SpecialFunctions.jl")
include("Strings.jl")
include("FileSave.jl")

simdb = nothing

#-------------------------------------------------------------------------------
"""
    function run(instance::SpeckleInstance)

Calculates the photon counts and correlations for a single instance of frequencies and phases.
Returns (SpeckleReadout,Correlation)
"""
function run(instance::SpeckleInstance)
    # calculate the max index of the time integration window
    nwindow = length(instance.γint)÷convert(Int, ceil(instance.params.duration/instance.params.window))

    readout = denseReadout(instance.γint,instance.bs,instance.detect)
    corr = correlate1d(readout,1,nwindow)
    return readout,corr
end

"""
    function run(params::SpeckleParams)

Calculates the photon counts and correlations for the number of repeats designated in params.
Returns a SpeckleSim object.
"""
function run(params::SpeckleParams)
    id = uuid4()
    dt = now()
    bs = Beamsplitter(1,1)
    readout = SpeckleReadout[]
    corr = Correlation[]
    sim = SpeckleSim(dt,id,params,bs,readout,corr)

    # make a directory to store results from this run
    results_data = joinpath(resultsDir(),string(id),"data")
    mkpath(results_data)
    results_plots = joinpath(resultsDir(),string(id),"plots")
    mkpath(results_plots)

    # instantiate frequencies and phases
    instance = SpeckleInstance(params)

    for i = 1:params.repeat
        # create photon counts and correlations
        readout, corr = run(instance)
        
        push!(sim.readout,readout)
        push!(sim.corr,corr)

        # reinstantiate frequencies and phases if desired
        if i != params.repeat && params.reinstance == true
            instance = SpeckleInstance(params)
        end
    end
    save(sim,results_data)
    return sim
end

function run(allparams::Dict; results_dir::String = results_directory)
    if results_dir != results_directory
        resultsDir(results_dir)
    end
    dbpath = joinpath(resultsDir(),"simdb.csv")
    if isfile(dbpath)
        global simdb = loadtable(dbpath)

    end
    paramVec = SpeckleParamsVector(allparams)
    simVec   = run.(paramVec)
    simTbl = tabulate(simVec)
    if simdb != nothing
        global simdb = mergeall(simdb,simTbl)
    else
        simdb = simTbl
    end
    JuliaDB.save(simdb,dbpath)
    return nothing
end
export run
#-------------------------------------------------------------------------------
end
