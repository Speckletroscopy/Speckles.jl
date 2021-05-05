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


"""
    fieldInstance!(efieldp::eFieldParams,n::Integer,bs::Beamsplitter,df::DataFrame;update::String = "cat")

Advances df by an instance of the intensity calculation
"""
function fieldInstance!(efieldp::eFieldParams,n::Integer,bs::Beamsplitter,df::DataFrame;update::String = "cat")
    eInstance = eFieldInstance(n,efieldp)
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

export fieldInstance!

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

export corrInstance!

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

export ftInstance!

"""
    classicalInstance!(efieldp::eFieldParams,n::Integer,bs::Beamsplitter,cuts::Tuple{T,T},idf::DataFrame,cdf::DataFrame,ftdf::DataFrame;update::String = "cat") where {T<:Integer}

Advances idf, cdf, and ftdf by instances of the intensity, correlation, and correlation Fourier transform respectively
"""
function classicalInstance!(efieldp::eFieldParams,n::Integer,bs::Beamsplitter,cuts::Tuple{T,T},idf::DataFrame,cdf::DataFrame,ftdf::DataFrame;update::String = "cat") where {T<:Integer}
    fieldInstance!(efieldp,n,bs,idf,update = update)
    corrInstance!(idf,cdf,update=update)
    ftInstance!(cdf,ftdf,cuts,update=update)
    return idf,cdf,ftdf
end

export classicalInstance!

end
