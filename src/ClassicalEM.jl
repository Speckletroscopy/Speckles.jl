################################################################################
# Electric field parameters
################################################################################
"""
    eFieldParams(Em::Vector,ωm::Vector,ω0::Number,σ::Number,rng::MersenneTwister)

Static parameters for the electric field
"""
struct eFieldParams
    Em::Vector # field magnitudes
    ωm::Vector # central frequencies of lines
    ω0::Number # central frequency for composite line
    σ::Number # Doppler broadening
    rng::MersenneTwister
    function eFieldParams(Em::Vector,ωm::Vector,ω0::Number,σ::Number,rng::MersenneTwister)
        @assert length(Em) == length(ωm) "Length of field magnitude and frequency vectors must match"
        new(
            convert(Vector{Complex},Em),
            convert(Vector{Real},ωm),
            convert(Real,ω0),
            convert(Real,σ),
            rng
        )
    end
end


"""
    eFieldParams(Em::Vector,ωm::Vector,σ::Number,seed::Integer=-1)

Static parameters for the electric field
"""
function eFieldParams(Em::Vector,ωm::Vector,σ::Number,seed::Integer = -1)
    if seed == -1
        rng = MersenneTwister()
    else
        rng = MersenneTwister(seed)
    end
    return eFieldParams(Em,ωm,ωm[1],σ,rng)
end

export eFieldParams

################################################################################
# Electric field instantiation
################################################################################
"""
	function seed!(params::eFieldParams,seed::Integer)

Reset the seed for the random number generator in params
"""
function seed!(params::eFieldParams,seed::Integer)
    if seed < 0
        Random.seed!(params.rng)
    else
        Random.seed!(params.rng,seed)
    end
end

export seed!

"""
Container holding frequencies and phases for one realization of the electric field.
"""
struct eFieldInstance
    ωn::Vector
    ϕmn::Matrix
    params::eFieldParams

    function eFieldInstance(ωn::Vector,ϕmn::Matrix,params::eFieldParams)
        @assert (size(ϕmn)[1] == length(params.ωm) && size(ϕmn)[2] == length(ωn)) "Phase array must have shape (m,n)"
        new(ωn,ϕmn,params)
    end
end

"""
    eFieldInstance(n::Integer,params::eFieldParams)

Generate frequencies and phases for the electric field from n emitting atoms and return an eFieldInstance object.
"""
function eFieldInstance(n::Integer,params::eFieldParams)
    ωDist = Normal(params.ω0,params.σ)
    ωn    = rand(params.rng,ωDist,n)
    ϕmn   = 2*π*rand(params.rng,Float64,(length(params.ωm),n))
    return eFieldInstance(ωn,ϕmn,params)
end


export eFieldInstance

"""
	function electricField(t::Real,params::eFieldParams)

Returns the electric field value at time t
"""
function electricField(t::Real,instance::eFieldInstance)
	Δm = instance.params.ωm .- instance.params.ω0
	# generate frequencies
	ωmn = transpose(instance.ωn) .+ Δm
	# add the phase
	exponentmn = -im*(t*ωmn+instance.ϕmn)
	# put them in the exponent
	enm = exp.(exponentmn)
	# multiply by the field magnitude
	fieldnm = instance.params.Em .* enm
	# add it all together
	return sum(ivec(fieldnm))
end

export electricField

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
    intensity(t::Real,instance::eFieldInstance)

Returns the intensity given an instance of the EM field and a time t.
"""
function intensity(t::Real,instance::eFieldInstance)
    eFieldT = electricField(t,instance)
    return intensity(eFieldT)
end

export intensity
