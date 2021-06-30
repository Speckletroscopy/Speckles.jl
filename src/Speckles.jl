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
@reexport using DSP
@reexport using Measurements
using CatViews

# location where all results are stored
results_directory = "results"

include("LightSource.jl")
include("Detector.jl")
include("Parameters.jl")
include("Correlation.jl")
include("SimulationFlow.jl")
include("Analysis.jl")
include("Strings.jl")
include("Plots.jl")
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
    @info "Beginning run id:$id"
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
    save(sim)
    return sim
end

"""
    function run(allparams::Dict; results_dir::String = results_directory)

Run a series of simulations and store them to results_dir.

- allparams is a dictionary with the same keys as SpeckleParams, but the
    values are lists with each simulation parameter to be run
"""
function run(allparams::Dict; results_dir::String = results_directory)
    # use the default name for the results directory if none is given
    if results_dir != results_directory
        resultsDir(results_dir)
    end

    # set the simulation database file path and check for existence
    dbpath = joinpath(resultsDir(),"simdb.csv")
    if isfile(dbpath)
        # load the simulation database file if it exists
        global simdb = JuliaDB.load(dbpath)
    else
        # set it to nothing if it doesn't exists so we can make a new one later
        global simdb = nothing
    end

    # split up allparams into a vector of SpeckleParams
    paramVec = SpeckleParamsVector(allparams)

    # run the simulation for each set of parameters
    simVec= run.(paramVec)

    # store the simulation results in a table
    simTbl = tabulate(simVec)

    # merge simulation results into database, if it exists
    #   make a new one if it doesn't...
    if simdb !== nothing
        global simdb = mergeall(simdb,simTbl)
    else
        simdb = simTbl
    end

    # *** ANALYSIS *** #
    # take fourier transform of correlations
    fftVec = SpeckleFFT.(simVec)

    # calculate snr
    snrVec = map(ft_par->snr(ft_par[1],ft_par[2]),zip(fftVec,paramVec))

    # function to summarize results of snr calculations in the simdb table
    function makeSnrDbEntry(t::IndexedTable,id::UUID)
        snrDict = Dict{Symbol,AbstractVector}()
        for col in colnames(t)
            if col == :type continue end
            if length(select(t,col))>1
                val = [mean(select(t,col))±std(select(t,col))]
            else
                val = select(t,col)
            end
            tname = select(t,:type)[1]
            nameStr = string("SNR_",tname,"Δ",col)
            nameSym = Symbol(nameStr)
            snrDict[nameSym] = val
        end
        snrDict[:id] = [id]
        out = table(snrDict,pkey=:id)
        @assert length(out) == 1 "Something went wrong"
        return out
    end

    allSnrEntries = nothing
    # TODO: make this more efficient
    for (i,snrTbl) in enumerate(snrVec)
        singleTbl = filter(t->t.type == "single",snrTbl)
        sumFFTtbl = filter(t->t.type == "sumFFT",snrTbl)
        FFTsumTbl = filter(t->t.type == "FFTsum",snrTbl)
        singleEntry = makeSnrDbEntry(singleTbl,simVec[i].id)
        sumFFTentry = makeSnrDbEntry(sumFFTtbl,simVec[i].id)
        FFTsumEntry = makeSnrDbEntry(FFTsumTbl,simVec[i].id)

        snrEntry = join(sumFFTentry,FFTsumEntry,how=:inner)
        snrEntry = join(singleEntry,snrEntry,how=:inner)
        if allSnrEntries === nothing
            allSnrEntries = snrEntry
        else
            allSnrEntries = mergeall(allSnrEntries,snrEntry)
        end
    end

    simdb = mergeall(simdb,allSnrEntries)

    # save the database to file and return the results
    JuliaDB.save(simdb,dbpath)
    return simdb
end
export run
#-------------------------------------------------------------------------------
end
