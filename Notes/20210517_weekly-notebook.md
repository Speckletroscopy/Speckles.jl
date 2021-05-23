#weekly #notes

| Parameter                | Value |
| ------------------------ | ----- |
| time resolution (ns)     | 0.01  |
| maximum time (ns)        | 20    |
| number of atoms          | 10    |
| Doppler broadening (GHz) | 10    |
| Number of lines          | 2     |                    |       |

## Classical Plots
### Instances: 1
![](../plots/archived/202105013_tres=0.01_tmax=10.0_len-%CE%BDM=2_ntot=100000_nbar=10_bigN=10_reset=1.0_temp=5000_len-mag=2_seed=-1_classical-single.svg)

### Instances: $10^4$


![](../plots/archived/202105013_tres=0.01_tmax=10.0_len-%CE%BDM=2_ntot=100000_nbar=10_bigN=10_reset=1.0_temp=5000_len-mag=2_seed=-1_classical-sum.svg)


## Photon Counting Plots

### Approach \#1
This approach treats the photon counting case similarly to the classical field case. The plot below was acquired by the following:
- Taking the intensity from one instance of frequencies and generating an average of 10 Poisson-distributed photons per instance in two beams
- Calculating the correlation between each beam 
$$
g_i = \sum_{j=0}^{N_t/2}\gamma_j\gamma_{i+j}	
$$

- where above, I define the following
	- $\gamma_j$ is the number of counts in bin $j$
	- $N_t$ is the number of time bins
- If $g_i$ is non-zero, I register one correlation count at a time of $t_{res}i$
- Repeat the above for each instance until the desired number of photons is reached
- Make a histogram of all time correlations
- Take the Fourier transform of the resulting histogram

![](../plots/archived/202105013_tres=0.01_tmax=10.0_len-νM=2_ntot=100000_nbar=10_bigN=10_reset=1.0_temp=5000_len-mag=2_seed=-1_photon-correlation.svg)

We can see that there is no visible correlation at the line-difference frequency. This is probably because (a) I use a boolean correlation function which wipes out the magnitude of each correlation instance, (b) the Poisson variance is large since the number of counts is small, and (c) the re-instancing the frequencies and phases does not allow the Poisson-distributed photon counts to approach an average.

### Approach \#2

We can see from the classical plots that the method in our paper *does* work in the classical limit. The problem is that we must gather enough photon counts to reduce the shot noise. However, due to detector deadtime we will be unable to acquire photon counts at the maximum time resolution continuously. Additionally, if we only correlate photon pairs within a single readout window, we get a nice lesson in the exponential nature of Poisson-distributed phenomena:

![](../plots/archived/202105013_tres=0.01_tmax=10.0_len-νM=2_ntot=100000_nbar=10_bigN=10_reset=1.0_temp=5000_len-mag=2_seed=-1_coincident-counts-vs-time.svg)

Obviously, we require simulation of detector deadtime to gather enough photon counts. This past week, I've implemented some new objects that will make it easier to simulate detector effects and light source properties:

```julia
"""  	Detector(deadtime::Number,resolution::Number,jitter::Number,efficiency::Number,darkcounts::Number)

Returns a Detector object.

Inputs:
    deadtime   : dead time in nanoseconds
    resolution : time resolution in nanoseconds
    jitter     : timing jitter in nanoseconds
    efficiency : quantum efficiency ∈ [0,1]
    darkcounts : dark count rate in GHz
"""
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
```


```julia
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
```

All instances of the EM field can be instantiated from a single LightSource object

```Julia
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
```

where the eField object is given by

```julia
"""
    eField(νn::Vector,ϕmn::Matrix,source::LightSource)

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
```

Photon counts are generated by first generating the average photon count rate over a series of time bins. This is accomplished by feeding the expected number of photon counts, nbar, a time duration t, a time resolution dt, and the above eField object to the following

```julia
function γIntensity(nbar::Real,t::Real,dt::Real,field::eField)
    intensitySeries = map(time->intensity(time,field),0:dt:t)
    return γIntensity(nbar,intensitySeries)
end
```

which returns a $\gamma$Intensity object

```julia
"""
    γIntensity(nbar::T,γvec::Vector{T}) where T<:Real

Returns γIntensity object which contains the average photon count rate and the renormalized intensity time series
"""
struct γIntensity{T<:Real}
    nbar::T
    γvec::Vector{T}
    function γIntensity(nbar::T,γvec::Vector{T}) where T<:Real
        total = sum(γvec)
        coeff = nbar/total
        return new{T}(nbar,coeff*γvec)
    end
end
```

Photon counts are then read out with the "readout" function

```julia
"""
    readout(γint:γIntensity,detect::Detector)

Returns a sparse vector containing the nonzero counts received by the detector for the given photon intensity time-series
"""
function readout(γint::γIntensity,detect::Detector)
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
```

which returns a sparse array of photon counts with the same time duration and resolution as the input $\gamma$Intensity object. Currently, deadtime is the only detector effect which is implemented. 

In development is an algorithm that will quickly calculate the correlation between two sparse readout vectors, one for each beam, which does not rely on sequentially iterating through values which will return zero.