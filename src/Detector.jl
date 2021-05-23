struct Detector
    deadtime::Float64
    resolution::Float64
    jitter::Float64
    efficiency::Float64
    darkcounts::Float64
    function Detector(deadtime::Float64,resolution::Float64,jitter::Float64,efficiency::Float64,darkcounts::Float64)
        @assert 0 ≤ efficiency ≤ 1 "Efficiency must be between 0 and 1"
        new(
            deadtime,
            resolution,
            jitter,
            efficiency,
            darkcounts
        )
    end
end

"""
    Detector(deadtime::Number,resolution::Number,jitter::Number,efficiency::Number,darkcounts::Number)

Returns a Detector object.

Inputs:
    deadtime   : dead time in nanoseconds
    resolution : time resolution in nanoseconds
    jitter     : timing jitter in nanoseconds
    efficiency : quantum efficiency ∈ [0,1]
    darkcounts : dark count rate in GHz
"""
function Detector(deadtime::Number,resolution::Number,jitter::Number,efficiency::Number,darkcounts::Number)
    return Detector(
        convert(Float64,deadtime),
        convert(Float64,resolution),
        convert(Float64,jitter),
        convert(Float64,efficiency),
        convert(Float64,darkcounts)
    )
end

export Detector

"""
    sparseReadout(t::Number,source::LightSource,detect::Detector)

Returns a sparse vector containing the nonzero counts received by the detector for the given duration t and LightSource
"""
function sparseReadout(t::Number,source::LightSource,detect::Detector)
    field = eField(source)
    return sparseReadout(t,field,detect)
end

"""
    sparseReadout(t::Number,field::eField,detect::Detector)

Returns a sparse vector containing the nonzero counts received by the detector for the given duration t and EM field
"""
function sparseReadout(t::Number,field::eField,detect::Detector)
    nb = nbar(t,field.source)
    γint = γIntensity(nb,t,detect.resolution,field)
    return sparseReadout(γint,detect)
end

"""
    sparseReadout(γint:γIntensity,detect::Detector)

Returns a sparse vector containing the nonzero counts received by the detector for the given photon intensity time-series
"""
function sparseReadout(γint::γIntensity,detect::Detector)
    out = sparsevec(Integer[],Integer[],length(γint))
    ideadtime = convert(Int,ceil(detect.deadtime/detect.resolution))
    i = 1
    while i ≤ length(γint)
        ct = poissonCount(γint[i])
        if ct != 0
            out[i] = ct
            i+=ideadtime
        end
        i+=1
    end
    return out
end

export sparseReadout

################################################################################
"""
    denseReadout(t::Number,source::LightSource,detect::Detector)

Returns a dense vector containing the counts received by the detector for the given duration t and LightSource
"""
function denseReadout(t::Number,source::LightSource,detect::Detector)
    field = eField(source)
    return denseReadout(t,field,detect)
end

"""
    denseReadout(t::Number,field::eField,detect::Detector)

Returns a dense vector containing the counts received by the detector for the given duration t and EM field
"""
function denseReadout(t::Number,field::eField,detect::Detector)
    nb = nbar(t,field.source)
    γint = γIntensity(nb,t,detect.resolution,field)
    return denseReadout(γint,detect)
end

"""
    denseReadout(γint:γIntensity,detect::Detector)

Returns a dense vector containing the counts received by the detector for the given photon intensity time-series
"""
function denseReadout(γint::γIntensity,detect::Detector)
    rngvec = Poisson.(γint.γvec)
    out = zeros(Int,length(γint))
    ideadtime = convert(Int,ceil(detect.deadtime/detect.resolution))
    i = 1
    while i ≤ length(γint)
        ct = rand(rngvec[i])
        if ct != 0
            out[i] = ct
            i+=ideadtime
        end
        i+=1
    end
    return out
end

export denseReadout