export plotsDir
export dataDir
export Beamsplitter
export classicalFreqData
export classicalSinglePlot
export classicalSumPlot

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
    γCorrTimePlot(timeVec::Vector, g2Vec::Vector, params::Dict, prefix::String)

Saves plot of photon correlation time series
"""
function γCorrTimePlot(timeDF::DataFrame, params::Dict, prefix::String)
    γplotName = string(prefix,"time-domain-photon-correlation.svg")
    inzcorr  = timeDF[!,:corr1] .> 0
    nzcorr = timeDF[!,:corr1][inzcorr]
    nztime = timeDF[!,:time][inzcorr]
	γplot = plot(nztime,nzcorr,label = false)
    xlabel!(L"\tau \textrm{ (ns)}")
	ylabel!("\$g^{(2)}(\\tau)\$")
    title!("Photon Correlations vs Time Offset")
    savefig(γplot,plotsDir(γplotName))
    @info "Saved $(plotsDir(γplotName))"
    if params[:repeat] > 1
        γplotName = string(prefix,"time-domain-photon-correlation-sum.svg")
        nzcorr = timeDF[!,:sum][inzcorr]
        γplot = plot(nztime,nzcorr,label = false)
        xlabel!(L"\tau \textrm{ (ns)}")
        ylabel!("\$\\sum_{i=1}^{$(params[:repeat])}g_i^{(2)}(\\tau)\$")
        title!("Sum of Photon Correlations vs Time Offset")
        savefig(γplot,plotsDir(γplotName))
        @info "Saved $(plotsDir(γplotName))"
    end
end

export γCorrTimePlot

"""
    γCorrFreqPlot(freqVec::Vector, fftVec::Vector, params::Dict, prefix::String)

Saves plot of photon correlation Fourier transform
"""
function γCorrFreqPlot(freqDF::DataFrame, params::Dict, prefix::String)
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
    savefig(γplot,plotsDir(γplotName))
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
        savefig(γplot,plotsDir(γplotName))
    end
end

export γCorrFreqPlot