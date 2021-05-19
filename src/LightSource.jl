
struct LightSource
    n::Integer # number of atoms
    Em::Vector # line magnitudes
    νm::Vector # central frequencies of lines
    σ::Number # Doppler broadening
    νMin::Number # minimum bandpass frequency
    νMax::Number # maximum bandpass frequency
    γRate::Number # photon count rate in Hertz
    function LightSource(n::Integer,Em::Vector,νm::Vector,σ::Number,νMin::Number,νMax::Number,γRate::Number)
        @assert length(Em) == length(νm) "Vectors for line magnitudes and frequencies must have the same length"
        @assert νMin < νMax "Bandpass minimum must be less than maximum"
        new(
            n,
            convert(Vector{Complex},Em),
            convert(Vector{Real},νm),
            convert(Real,σ),
            convert(Real,νMin),
            convert(Real,νMax),
            convert(Real,γRate)
        )
    end
end

function LightSource(n::Integer,Em::Vector,νm::Vector,σ::Number, γRate::Number)
    νMin = min(νm) - 5*σ # automatically set the bandpass minimum to 5σ below the lowest line frequency
    νMax = max(νm) + 5*σ # automatically set the bandpass maximum to 5σ above the highest line frequency
    return LightSource(n,Em,νm,σ,νMin,νMax,γRate)
end

export LightSource

"""
    ν0(source::LightSource)

Returns the average frequency of spectral lines in source
"""
function ν0(source::LightSource)
    return mean(source.νm)
end

export ν0

"""
    Δm(source::LightSource)

Returns vector of the sources line separation from the average frequency of all lines
"""
function Δm(source::LightSource)
    return source.νm .- ν0(source) 
end

export Δm

"""
    lineShifts(source::LightSource)

Returns a vector with the frequency differences between all spectral lines
"""
function lineShifts(source::LightSource)
    return map(x->abs(x[2]-x[1]),subsets(source.νm,2))
end

export lineShifts

################################################################################
# Electric field 
################################################################################
"""
Container holding frequencies and phases for one realization of the electric field.
"""
struct eField
    νn::Vector
    ϕmn::Matrix
    source::LightSource

    function eField(ωn::Vector,ϕmn::Matrix,source::LightSource)
        @assert (size(ϕmn)[1] == length(source.ωm) && size(ϕmn)[2] == length(ωn)) "Phase array must have shape (m,n)"
        new(ωn,ϕmn,source)
    end
end

"""
    eField(source::LightSource, seed::Integer = -1)

Generate frequencies and phases for a single instance of the electric field.
"""
function eField(source::LightSource, seed::Integer = -1)
    if seed != -1
       Random.seed!(seed) 
    end
    νDist = Normal(ν0(source),source.σ)
    νn    = rand(νDist,source.n)
    ϕmn   = 2*π*rand(Float64,(length(source.νm),source.n))
    return eField(νn,ϕmn,params)
end

export eField

"""
	function eField(t::Real,field::eField)

Returns the electric field value at time t
"""
function eFieldT(t::Real,field::eField)
	# generate frequencies
	ωmn = transpose(field.νn) .+ Δm(field.source)
    ωmn .*= 2*π
	# add the phase
	exponentmn = -im*(t*ωmn+field.ϕmn)
	# put them in the exponent
	enm = exp.(exponentmn)
	# multiply by the field magnitude
	fieldnm = field.source.Em .* enm
	# add it all together
	return sum(ivec(fieldnm))
end

export eFieldT

################################################################################
# Calculate intensity
################################################################################
"""
    intensity(eFieldT::Number)

Returns the intensity given the value of the electric field at a particular point in time and space
"""
function intensity(eFieldT::Number)
    return real(eFieldT*conj(eFieldT))
end

"""
    intensity(t::Real,field::eField)

Returns the intensity given an instance of the EM field and a time t.
"""
function intensity(t::Real,field::eField)
    eFieldT = eFieldT(t,field)
    return intensity(eFieldT)
end

export intensity

"""
    Beamsplitter(r::Number,t::Number)

Struct to hold beam splitter coefficients.
"""
struct Beamsplitter
    r::Number
    t::Number
    function Beamsplitter(r::Number,t::Number)
        beamNorm = sqrt(r^2+t^2)
        new(r/beamNorm,t/beamNorm)
    end
end


