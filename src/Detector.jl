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
    Detector(deadtime::Number,resolution::Number,efficiency::Number,darkcounts::Number)

Returns a Detector object.

Inputs:
    deadtime   : dead time in nanoseconds
    resolution : time resolution in nanoseconds
    jitter     : timing jitter in nanoseconds
    efficiency : quantum efficiency ∈ [0,1]
    darkcounts : dark count rate
"""
function Detector(deadtime::Number,resolution::Number,efficiency::Number,darkcounts::Number)
    return Detector(
        convert(Float64,deadtime),
        convert(Float64,resolution),
        convert(Float64,jitter),
        convert(Float64,efficiency),
        convert(Float64,darkcounts)
    )
end

"""
    readout(t::Number,source::LightSource,detect::Detector)

Returns the counts received by a detector from the given source over a time t (ns)
"""
function readout(t::Number,source::LightSource,detect::Detector)
    field = eField(source)
    return readout(t,field,detect)
end

"""
    readout(t::Number,field::eField,detect::Detector)

Returns the counts received by a detector from the given incident EM field over a time t (ns)
"""
function readout(t::Number,field::eField,detect::Detector)
    nbar = t*field.source.γRate # mean number of photons expected over time t
    intensityVec = map(time->intensity(time,field),0:detect.resolution:t)
    γint = γIntensity(nbar,intensityVec)
    return readout(γint,detect)
end

function readout(γint:γIntensity,detect::Detector)

end

function generate(t::Number,dt::Number,dead::Integer,γRate::Number)
    nbar = γRate*dt
    out = Integer[]
    i = 1
    while (i-1)*dt < t
        if poissonCount(nbar) != 0
            push!(out,i)
            i+=dead
        else
            i+=1
        end
    end
    return out
end

export generate
