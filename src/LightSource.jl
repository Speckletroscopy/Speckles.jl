
struct LightSource
    n::Integer # number of atoms
    Em::Vector # line magnitudes
    νm::Vector # central frequencies of lines in GHz
    σ::Number # Doppler broadening in GHz
    νMin::Number # minimum bandpass frequency in GHz
    νMax::Number # maximum bandpass frequency in GHz
    γRate::Number # photon count rate in GHz
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
    νMin = min(νm...) - 5*σ # automatically set the bandpass minimum to 5σ below the lowest line frequency
    νMax = max(νm...) + 5*σ # automatically set the bandpass maximum to 5σ above the highest line frequency
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

    function eField(νn::Vector,ϕmn::Matrix,source::LightSource)
        @assert (size(ϕmn)[1] == length(source.νm) && size(ϕmn)[2] == length(νn)) "Phase array must have shape (m,n)"
        new(νn,ϕmn,source)
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
    return eField(νn,ϕmn,source)
end

export eField

"""
	function eField(t::Number,field::eField)

Returns the electric field value at time t(ns)
"""
function eFieldT(t::Number,field::eField)
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
    intensity(et::Number)

Returns the intensity given the value of the electric field at a particular point in time and space
"""
function intensity(et::Number)
    return real(et*conj(et))
end

"""
    intensity(t::Real,field::eField)

Returns the intensity given an instance of the EM field and a time t(ns).
"""
function intensity(t::Number,field::eField)
    efield = eFieldT(t,field)
    return intensity(efield)
end

export intensity

"""
    γIntensity{T}(nbar::T,intensity::Vector{T}) where {T<:Real}

Returns γIntensity object which contains the average photon count rate and the renormalized intensity time series
"""
struct γIntensity{T<:Real}
    nbar::T
    intensity::Vector{T}
    function γIntensity{T}(nbar::T,intensity::Vector{T}) where {T<:Real}
        total = sum(intensity)
        coeff = nbar/total
        nintensity = intensity .* coeff
        new(nbar,nintensity)
    end
end

export γIntensity

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


