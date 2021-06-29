################################################################################
# File tree manipulation
################################################################################
#-------------------------------------------------------------------------------
"""
    resultsDir()

Returns the name of the directory where results are being stored
"""
function resultsDir()
    return results_directory
end

"""
    resultsDir(newdir::String)

Sets the directory where results are stored to newdir
"""
function resultsDir(newdir::String)
    global results_directory = newdir
    return results_directory
end
export resultsDir

#-------------------------------------------------------------------------------

"""
    plotsDir(sim::SpeckleSim)

Returns the path to the plots directory for the input simulation
"""
function plotsDir(sim::SpeckleSim)
    return joinpath(resultsDir(),string(sim.id),"plots")
end
export plotsDir

#-------------------------------------------------------------------------------

"""
    dataDir(name::String,dirname::String = "data")

Concatenates string with name of data directory
"""
function dataDir(sim::SpeckleSim)
    return joinpath(resultsDir(),string(sim.id),"data")
end
export dataDir

#-------------------------------------------------------------------------------
################################################################################
# Database management
################################################################################
"""
    mergeall(a::IndexedTable,b::IndexedTable)

Merges tables a and b. Unmatched columns are filled in with missing values.
"""
function mergeall(a::IndexedTable,b::IndexedTable)
    # make column names into sets for set operations
    colseta = Set(colnames(a))
    colsetb = Set(colnames(b))
    allcols = union(colseta, colsetb)
    commoncols = intersect(colseta,colsetb)
    # return just the regular merge if all columns are shared
    if commoncols == allcols
        return merge(a,b)
    end

    # get column names not common between the tables
    cols_anotb = setdiff(colseta,colsetb)
    cols_bnota = setdiff(colsetb,colseta)

    # merge parts of table with common columns
    commona = select(a,Tuple(commoncols))
    commonb = select(b,Tuple(commoncols))
    commontbl = merge(commona,commonb)
    # join uncommon columns into the table
    if length(cols_anotb) > 0
        onlya = select(a,(:id,cols_anotb...))
        commontbl = join(commontbl,onlya; how=:outer)
    end
    if length(cols_bnota) > 0
        onlyb = select(b,(:id,cols_bnota...))
        commontbl = join(commontbl,onlyb; how=:outer)
    end
    return commontbl
end

#-------------------------------------------------------------------------------
function tabulate(sim::SpeckleSim)
    simDict = Dict{Symbol,AbstractVector}()
    # iterate through fields in sim.params and expand arrays to individual columns
    for key in fieldnames(typeof(sim.params))
        if typeof(getfield(sim.params,key)) <: AbstractArray
            keyStr = string(key)
            for (i,val) in enumerate(getfield(sim.params,key))
                iKeyStr = string(keyStr,i)
                simDict[Symbol(iKeyStr)] = [val]
            end
        else
            simDict[key] = [getfield(sim.params,key)]
        end
    end
    simDict[:dt]  = [sim.dt]
    simDict[:id]  = [sim.id]
    simDict[:bst] = [sim.bs.t]
    simDict[:bsr] = [sim.bs.r]
    return table(simDict; pkey=:id)
end

function tabulate(simvec::Vector{T}) where {T<:SpeckleSim}
    out = tabulate(simvec[1])
    if length(simvec) == 1
        return out
    else
        for i=2:length(simvec)
            simtbl = tabulate(simvec[i])
            out = mergeall(out,simtbl)
        end
    end
    return out
end
export tabulate
#-------------------------------------------------------------------------------

################################################################################
# Save functions
################################################################################
function save(sim::SpeckleSim)
    for i = 1:length(sim.readout)
        
        beamDict = Dict(:b1=>sim.readout[i].beam1, :b2=>sim.readout[i].beam2)
        beamTbl  = table(beamDict)
        beamName = string("counts",i,".csv")
        beamPath = joinpath(dataDir(sim),beamName)
        JuliaDB.save(beamTbl,beamPath)

        corrDict = Dict{Symbol,Vector}()
        if typeof(sim.corr[i]) <: CorrelationVector
            corrDict[Symbol(sim.corr[i].n)] = sim.corr[i].data
        else
            for j=1:size(sim.corr[i].data)[2]
                corrDict[Symbol(j)] = sim.corr[i].data[:,j]
            end
        end
        corrTbl = table(corrDict)
        corrName = string("correlation",i,".csv")
        corrPath = joinpath(dataDir(sim),corrName)
        JuliaDB.save(corrTbl,corrPath)
    end
    @info "Saved data for simulation id: $(sim.id)"
    return nothing
end

#-------------------------------------------------------------------------------
function save(plt::Plots.Plot,name::String,sim::SpeckleSim)
    figPath = joinpath(plotsDir(sim),name)
    savefig(plt,figPath)
    @info "Saved $name in $(plotsDir(sim))"
end
#-------------------------------------------------------------------------------
"""
    γIntensityPlot(times::Array,γint::Array,prefix::String)

Plots average photon counts and saves to file. Returns relative path to resulting plot.
"""
function γIntensityPlot(times::Array,γint::Array,prefix::String)
    γplotName = string(prefix,"avg-photons-vs-time.svg")
    γplot = plot(times,γint,label=false)
    xlabel!("time (ns)")
    ylabel!("Avg photon counts")
    title!("Average Photon Counts vs Time")
    savefig(γplot,γplotName)
    @info "Saved $γplotName"
    return γplotName
end

export γIntensityPlot
#-------------------------------------------------------------------------------

"""
    γCountPlot(times::Array,γint::Array,prefix::String)

Plots 'measured' photon counts and saves to file. Returns relative path to resulting plot.
"""
function γCountPlot(times::Array,γcount::Array, γint::Array,prefix::String)
    γplotName = string(prefix,"photon-counts-vs-time.svg")
    γplot = bar(times,γcount,label="Simulated Counts")
    sig1 = sqrt.(γint)
    sig1low = map(x -> x < 1 ? x : sqrt(x),γint)
    plot!(γplot,times,γint,label = "Average Counts"; ribbon = (sig1low,sig1))
    xlabel!("time (ns)")
    ylabel!("photon counts")
    title!("Photon Counts vs Time")
    savefig(γplot,γplotName)
    @info "Saved $γplotName"
    return γplotName
end

export γCountPlot
"""
    γCorrTimePlot(timeVec::Vector, g2Vec::Vector, params::Dict, prefix::String)

Saves plot of photon correlation time series
"""
function γCorrTimePlot(timeDF::IndexedTable, params::Dict, prefix::String)
    out = []
    γplotName = string(prefix,"time-domain-photon-correlation.svg")
    inzcorr  = timeDF[!,:corr1] .> 0
    nzcorr = timeDF[!,:corr1][inzcorr]
    nztime = timeDF[!,:time][inzcorr]
	γplot = plot(nztime,nzcorr,label = false)
    xlabel!(L"\tau \textrm{ (ns)}")
	ylabel!("\$g^{(2)}(\\tau)\$")
    title!("Photon Correlations vs Time Offset")
    savefig(γplot,γplotName)
    push!(out,γplotName)
    @info "Saved $(γplotName)"
    if params[:repeat] > 1
        γplotName = string(prefix,"time-domain-photon-correlation-sum.svg")
        nzcorr = timeDF[!,:sum][inzcorr]
        γplot = plot(nztime,nzcorr,label = false)
        xlabel!(L"\tau \textrm{ (ns)}")
        ylabel!("\$\\sum_{i=1}^{$(params[:repeat])}g_i^{(2)}(\\tau)\$")
        title!("Sum of Photon Correlations vs Time Offset")
        savefig(γplot,γplotName)
        push!(out,γplotName)
        @info "Saved $(γplotName)"
    end
    return out
end

export γCorrTimePlot
#-------------------------------------------------------------------------------

"""
    γCorrFreqPlot(freqVec::Vector, fftVec::Vector, params::Dict, prefix::String)

Saves plot of photon correlation Fourier transform
"""
function γCorrFreqPlot(freqDF::IndexedTable, params::Dict, prefix::String)
    out = []
    γplotName = string(prefix,"frequency-domain-photon-correlation.svg")
	γplot = plot(freqDF[!,:freq],freqDF[!,:corr1],label = false)
	xlabel!("frequency (GHz)")
	ylabel!("\$\\hat{g}^{(2)}(\\nu)\$")
	title!("Photon Correlations vs Frequency")
    if length(params[:νm]) > 1
        shifts = map(x->abs(x[2]-x[1]),subsets(params[:νm],2))
        shiftname = length(params[:νm]) > 2 ? "frequency shifts" : "frequency shift"
        vline!(shifts,label=shiftname,ls=:dash)
    end
    savefig(γplot,γplotName)
    push!(out,γplotName)
    @info "Saved $(γplotName)"
    if params[:repeat] > 1
        γplotName = string(prefix,"frequency-domain-photon-correlation-sum.svg")
        γplot = plot(freqDF[!,:freq],freqDF[!,:sum],label = false)
        xlabel!("frequency (GHz)")
        ylabel!("\$\\sum_{i=1}^{$(params[:repeat])}\\hat{g}^{(2)}(\\nu)\$")
        title!("Sum of Photon Correlations vs Frequency")
        if length(params[:νm]) > 1
            shifts = map(x->abs(x[2]-x[1]),subsets(params[:νm],2))
            shiftname = length(params[:νm]) > 2 ? "frequency shifts" : "frequency shift"
            vline!(shifts,label=shiftname,ls=:dash)
        end
        savefig(γplot,γplotName)
        push!(out,γplotName)
        @info "Saved $(γplotName)"
    end
    return out
end

export γCorrFreqPlot
#-------------------------------------------------------------------------------
