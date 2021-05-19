export plotsDir
export dataDir
export Beamsplitter
export classicalFreqData
export classicalSinglePlot
export classicalSumPlot
export γCorrFreqPlot

"""
    plotsDir(name::String,dirname::String = "plots")

Concatenates string with name of plots directory
"""
function plotsDir(name::String,dirname::String = "plots")
    return joinpath(dirname,name)
end

"""
    dataDir(name::String,dirname::String = "data")

Concatenates string with name of data directory
"""
function dataDir(name::String,dirname::String = "data")
    return joinpath(dirname,name)    
end
"""
    classicalFreqData(dfFreqs::DataFrame,params::Dict,prefix::String)

Saves frequency dataframe to a csv file.
"""
function classicalFreqData(dfFreqs::DataFrame,params::Dict,prefix::String)
    dataName = string(prefix,"classical-frequency.csv")
    CSV.write(dataDir(dataName),dfFreqs)
end

"""
    classicalSinglePlot(dfFreqs::DataFrame, params::Dict, prefix::String)

Saves plot of single instance classical correlation Fourier transform.
"""
function classicalSinglePlot(dfFreqs::DataFrame, params::Dict, prefix::String)
    shifts = map(x->abs(x[2]-x[1]),subsets(params[:νM],2))
    singleplotName = string(prefix,"classical-single.svg")
    single = dfFreqs.PowerSpec
    ftSingle, ftFreqs = fftPositiveFreq(single,dfFreqs.freq)
    singleplot = plot(ftFreqs,abs.(ftSingle),label=false)
    vline!(shifts,label="Frequency shift",ls=:dash)
    xlabel!("frequency (GHz)")
	title!("\$\\hat{g}^{(2)}(\\nu)\$")
    savefig(singleplot,plotsDir(singleplotName))
    return nothing
end

"""
    classicalSumPlot(dfFreqs::DataFrame, params::Dict, prefix::String)

Saves plot of multi-look classical correlation Fourier transform.
"""
function classicalSumPlot(dfFreqs::DataFrame, params::Dict, prefix::String)
    shifts = map(x->abs(x[2]-x[1]),subsets(params[:νM],2))
    nInstances = convert(Integer,ceil(params[:ntot]/params[:nbar]))
    sumplotName = string(prefix,"classical-sum.svg")
    ftSum, ftFreqs = fftPositiveFreq(dfFreqs.sum,dfFreqs.freq)
    sumplot = plot(ftFreqs,abs.(ftSum),label=false)
    vline!(shifts,label="Frequency shift",ls=:dash)
    xlabel!("frequency (GHz)")
    title!("\$\\sum_{i=1}^{$(nInstances)}\\hat{g}_i^{(2)}(\\nu)\$")
    savefig(sumplot,plotsDir(sumplotName))
    return nothing
end

"""
    γCorrFreqPlot(freqVec::Vector, fftVec::Vector, params::Dict, prefix::String)

Saves plot of photon correlation Fourier transform
"""
function γCorrFreqPlot(freqVec::Vector, fftVec::Vector, params::Dict, prefix::String)
    shifts = map(x->abs(x[2]-x[1]),subsets(params[:νM],2))
    γplotName = string(prefix,"photon-correlation.svg")
	γplot = plot(freqVec,abs.(fftVec),label = false)
	xlabel!("frequency (GHz)")
	title!("Fourier transform of photon correlations")
	vline!(shifts,label="Frequency shift",ls=:dash)
    savefig(γplot,plotsDir(γplotName))
end
