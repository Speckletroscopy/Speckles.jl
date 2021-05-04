"""
    fftPositiveFreq(ft,f)

Returns (coefficients, frequencies) for positive frequencies of the fourier transform
"""
function fftPositiveFreq(fft,f)
	@assert length(fft) == length(f) "Fourier transform and frequency vectors much have matching length"
	fftRange = f .>= 0
	return (fft[fftRange],f[fftRange])
end

export fftPositiveFreq

"""
    meanFFT(fft::Vector{S},cuts::Tuple{T,T}) where {S<:Real, T<:Integer}

Takes mean before taking fourier transform
"""
function meanFFT(timeSeries::Vector{S},cuts::Tuple{T,T}) where {S<:Real, T<:Integer}
	tsPrep = timeSeries[1+cuts[1]:end-cuts[2]]
	tsPrepMean = mean(tsPrep)
	tsPrep = tsPrep .- tsPrepMean
	return fft(tsPrep)
end

export meanFFT

"""
    fftFreq(tres::Real,τ::Vector,cuts::Tuple{T,T}) where {T<:Integer}

Returns frequencies associated with Fourier transform 
"""
function fftFreq(tres::Real,τ::Vector,cuts::Tuple{T,T}) where {T<:Integer}
	τlen = length(τ) - sum(cuts)
	freqFFT = fftfreq(τlen,1/tres)
end

export fftFreq
