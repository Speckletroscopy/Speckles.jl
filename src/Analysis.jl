################################################################################
# Fourier Transforms
################################################################################
"""
    fftPositiveFreq(f,ft)

Returns (frequencies, coefficients) for positive frequencies of the fourier transform
"""
function fftPositiveFreq(ft::AbstractVector{T},f::AbstractVector{T}) where {T<:Real}
	@assert length(ft) == length(f) "Fourier transform and frequency vectors much have matching length"
	ftRange = f .>= 0
	return (f[ftRange],ft[ftRange])
end

export fftPositiveFreq

"""
    meanFFT(fft::Vector{S},cuts::Tuple{T,T}) where {S<:Real, T<:Integer}

Takes mean before taking fourier transform
"""
function meanFFT(timeSeries::Vector{S},cuts::Tuple{T,T}) where {S<:Real, T<:Integer}
	tsPrep = timeSeries[cuts[1]:cuts[2]]
	tsPrepMean = mean(tsPrep)
	tsPrep = tsPrep .- tsPrepMean
	return FFTW.fft(tsPrep)
end

export meanFFT

"""
    fftFreq(tres::Real,τ::Vector,cuts::Tuple{T,T}) where {T<:Integer}

Returns frequencies associated with Fourier transform 
"""
function fftFreq(tres::Real,cuts::Tuple{T,T}) where {T<:Integer}
	τlen = cuts[2]-cuts[1]+1
	freqFFT = FFTW.fftfreq(τlen,1/tres)
end

export fftFreq

"""
    niceFFT(timeSeries::Vector{S},cuts::Tuple{T,T},params::SpeckleParams) where {S<:Real, T<:Integer}

Performs some light data processing for a nice looking correlation fft.
    - Does a low cut on the time series to get past the non-periodic part of g2τ
    - Subtracts the mean of the time series to minimize the constant term coefficient

Returns (frequencies, coefficients)
"""
function niceFFT(timeSeries::Vector{S},params::SpeckleParams) where {S<:Real}
    # low cut 5 efolds beyond the non-periodic part of g2τ
    τlow = sqrt(5)/params.σ
    ilow = convert(Int,ceil(τlow/params.resolution))
    cuts = (ilow,length(timeSeries))
    mfft = meanFFT(timeSeries,cuts)
    mfft = abs.(mfft)
    mfft = convert(Vector{Float64},mfft)
    freqs = fftFreq(params.resolution,cuts)
    return fftPositiveFreq(freqs,mfft)
end
export niceFFT

struct SpeckleFFT{V<:Vector}
    freqs::Vector
    singles::Vector{V}
    sumFFT::Vector
    FFTsum::Vector
end
function SpeckleFFT(sim::SpeckleSim{U,V}) where {U<:SpeckleReadout,V<:Correlation}
    corrFFTvec = map(corr->niceFFT(corr.data,sim.params)[2],sim.corr)
    sumCorrFFT = +(corrFFTvec...)
    corrSum = +(getproperty.(sim.corr,:data)...)
    freqs,corrSumFFT = niceFFT(corrSum,sim.params)

    return SpeckleFFT(freqs,corrFFTvec,corrSumFFT,sumCorrFFT)
end
export SpeckleFFT
################################################################################
# Linear Prediction
################################################################################
function LPCoeff(corr::CorrelationVector)
    p = 20
    return DSP.LPC.lpc(corr.data,p)
end
################################################################################
# Analysis functions
################################################################################

function snr(ft::Vector,shifts::Vector,ishifts::Vector)
            # slice out the peaks and put the noise in a single view
            noiseviews = map(slice->view(ft,UnitRange(slice[1],slice[2])),ishifts)
            noiseviews = CatView(noiseviews...)
            # noise statistics
            μ = mean(noiseviews)
            σ = std(noiseviews)
            # select the max value in each peak window as the peak height
            peaks = Float64[]
            for k=1:length(shifts)
                peakview = view(ft,UnitRange(ishifts[k][2],ishifts[k+1][1]))
                peakval = max(peakview...)
                push!(peaks,peakval)
            end
            # subtract the mean value of the noise from the peak height
            peaks = peaks .- μ
            # calculate snr
            peaks = peaks/σ
            # zip together shifts and snrs for reference
            return collect(zip(shifts,peaks))
end

"""
    snr(sfft::SpeckleFFT,params::SpeckleParams)

Calculates the signal-to-noise ratio for each peak in sfft
"""
function snr(sfft::SpeckleFFT,params::SpeckleParams)

    if length(params.νm) > 1
        halfwindow = 2
        # find frequency shifts
        shifts = map(x->abs(x[2]-x[1]),subsets(params.νm,2))

        # find fft indices where frequency shifts are supposed to be
        shifts = sort(shifts)
        # shiftSym = Symbol.(shifts)
        ishifts = Tuple{Int64,Int64}[]
        i = 1
        istart = i
        iend = i
        for shift in shifts
            # iterate until we hit a shift or run out of frequencies
            while i > length(sfft.freqs) || sfft.freqs[i] < shift 
                i+=1
            end

            # back up one index plus half the size of the cut around the peak
            itest = i-1-halfwindow

            # check if itest is valid
            if i > length(sfft.freqs)
                iend = length(sfft.freqs)
            elseif itest < istart
                iend = istart
            else
                iend = itest
            end
            # store range for later use
            push!(ishifts,(istart,iend))
            # new starting point is on the other side of the peak
            istart = iend + 2*halfwindow
        end
        @assert i < length(sfft.freqs) "Counting should have stopped before the final shift"
        # store final slice 
        push!(ishifts,(istart,length(sfft.freqs)))
        @assert length(ishifts) == length(shifts)+1 "Slicing went wrong!"
        
        # do the singles first
        singleSnr = map(single->snr(single,shifts,ishifts),sfft.singles)
        sumFFTsnr = snr(sfft.sumFFT,shifts,ishifts)
        FFTsumSnr = snr(sfft.FFTsum,shifts,ishifts)

        return singleSnr,sumFFTsnr,FFTsumSnr

    else
        return nothing
    end
end
