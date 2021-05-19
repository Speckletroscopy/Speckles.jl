module Speckles

using Reexport
@reexport using Plots
@reexport using Dates
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
include("LightSource.jl")
include("SpeckleFunctions.jl")
include("SpeckleFT.jl")
include("FileSave.jl")

export fieldInstance!
export corrInstance!
export ftInstance!
export classicalInstance!
export paramVector
export makeName


"""
    fieldInstance!(source::LightSource,bs::Beamsplitter,df::DataFrame;update::String = "cat")

Advances df by an instance of the intensity calculation
"""
function fieldInstance!(source::LightSource,bs::Beamsplitter,df::DataFrame;update::String = "cat")
    eInstance = eField(source)
    efieldt = map(t->electricField(t,eInstance),df.time)
    efieldtBeam = bs.t * efieldt
    inty = intensity.(efieldtBeam)
    if update == "sum"
        df[!,"intensity"] = inty
if "sum" in names(df)
            df[!,"sum"] .+= inty
        else
            df[!,"sum"] = inty
        end
    else
        iname = string("intensity",size(df)[2])
        df[!,iname] = inty
    end
    return df
end

"""
    corrInstance!(idf::DataFrame,cdf::DataFrame;update::String = "cat")

Advances cdf by an instance of the corrlation calculation.
"""
function corrInstance!(idf::DataFrame,cdf::DataFrame;update::String = "cat")
    iname = string("g2tau",size(cdf)[2])
    if update == "sum"
        intensityBeam = idf.intensity
    else
        intensityBeam = idf[!,end]
    end
    window = size(cdf)[1]
    g2τNorm = mean(intensityBeam)^2
    g2τInst = map(0:window-1) do i
        autocorrelate(intensityBeam,i,window)/g2τNorm
    end
    if update == "sum"
        cdf[!,"g2tau"] = g2τInst
        if "sum" in names(cdf)
            cdf[!,"sum"] .+= g2τInst
        else
            cdf[!,"sum"] = g2τInst
        end
    else
        iname = string("g2tau",size(cdf)[2])
        cdf[!,iname] = g2τInst
    end
    return cdf
end

"""
    ftInstance!(cdf::DataFrame,ftdf::DataFrame,cuts::Tuple{T,T};update::String = "cat") where {T<:Integer}

Advances ftdf by an instance of the Fourier transform of the correlation calculation
"""
function ftInstance!(cdf::DataFrame,ftdf::DataFrame,cuts::Tuple{T,T};update::String = "cat") where {T<:Integer}
    if update == "sum"
        fft = meanFFT(cdf.g2tau,cuts)
        ftdf[!,"PowerSpec"] = fft
        if "sum" in names(ftdf)
            ftdf[!,"sum"] .+= fft
        else
            ftdf[!,"sum"] = fft
        end
    else
        fft = meanFFT(cdf[!,end],cuts)
        iname = string("PowerSpec",size(ftdf)[2])
        ftdf[!,iname] = fft
    end
    return ftdf
end

"""
    classicalInstance!(source::LightSource,n::Integer,bs::Beamsplitter,cuts::Tuple{T,T},idf::DataFrame,cdf::DataFrame,ftdf::DataFrame;update::String = "cat") where {T<:Integer}

Advances idf, cdf, and ftdf by instances of the intensity, correlation, and correlation Fourier transform respectively
"""
function classicalInstance!(source::LightSource,bs::Beamsplitter,cuts::Tuple{T,T},idf::DataFrame,cdf::DataFrame,ftdf::DataFrame;update::String = "cat") where {T<:Integer}
    fieldInstance!(source,bs,idf,update = update)
    corrInstance!(idf,cdf,update=update)
    ftInstance!(cdf,ftdf,cuts,update=update)
    return idf,cdf,ftdf
end

# function filenameGenerator
# function γCorrelate(intensity::Vector,nbar::Number,deadtime::Bool)
    # # calculate the average photon count rate per bin in each beam
    # nbarBeam = nbar/2
    # γAvgBeam = γIntensity(intensity,nbarBeam)

    # # generate counts for each beam
    # γCountsBeam1 = poissonCount.(γAvgBeam)
    # γCountsBeam2 = poissonCount.(γAvgBeam)

    # # 
# end

"""
    paramVector(params::Dict)

Splits multiple parameter dictionary into a vector of individual parameter dictionaries.
"""
function paramVector(params::Dict)
    k = collect(keys(params))
    v = collect(Iterators.product(collect(values(params))...))
    return map(vals->Dict(collect(zip(k,vals))),ivec(v))
end

"""
    makeName(params::Dict; excpt::Dict = Dict())

Returns a string based on parameters in the simulation dictionary.
"""
function makeName(params::Dict; excpt::Dict = Dict())
    date = today()
    mm = length(month(date)) == 1 ? string(0,month(date)) : month(date)
    dd = length(day(date)) == 1 ? string(0,day(date)) : day(date)
    out = string(year(date),mm,dd,"_")
    for key in keys(params)
        keyname = key
        val = params[key]
        if isa(val,Array)
            val = length(val)
            keyname = string("len-",keyname)
        end
        out = string(out,keyname,"=",val,"_")
    end
    return out
end

end