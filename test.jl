module Tst

using Speckles
using Plots


function fieldInstance!(efieldp::eFieldParams,n::Integer,bs::Beamsplitter,df::DataFrame;update::String = "cat")
    eInstance = eFieldInstance(n,efieldp)
    efieldt = map(t->electricField(t,eInstance),df.time)
    efieldtBeam = bs.t * efieldt
    inty = intensity.(efieldtBeam)
    if update == "sum"
        df[!,"intensity"] = inty
        if "sum" in names(df)
            df[!,"sum"] += inty
        else
            df[!,"sum"] = inty
        end
    else
        iname = string("intensity",size(df)[2])
        df[!,iname] = inty
    end
    return df
end


function corrInstance!(idf::DataFrame,cdf::DataFrame;update::String = "cat")
    iname = string("g2tau",size(cdf)[2])
    intensityBeam = idf[!,end]
    window = size(cdf)[1]
    g2τNorm = mean(intensityBeam)^2
    g2τInst = map(0:window-1) do i
        autocorrelate(intensityBeam,i,window)/g2τNorm
    end
    if update == "sum"
        if "sum" in names(cdf)
            cdf[!,"sum"] += g2τInst
        else
            cdf[!,"sum"] = g2τInst
        end
        cdf[!,"g2tau"] = g2τInst
    else
        iname = string("g2tau",size(cdf)[2])
        cdf[!,iname] = g2τInst
    end
    return cdf
end


function ftInstance!(cdf::DataFrame,ftdf::DataFrame,cuts::Tuple{T,T};update::String = "cat") where {T<:Integer}
    fft = meanFFT(cdf[!,end],cuts)
    if update == "sum"
        if "sum" in names(ftdf)
            ftdf[!,"sum"] += fft
        else
            ftdf[!,"sum"] = fft
        end
        ftdf[!,"PowerSpec"] = fft
    else
        iname = string("PowerSpec",size(ftdf)[2])
        ftdf[!,iname] = fft
    end
    return ftdf
end

function classicalInstance!(efieldp::eFieldParams,n::Integer,bs::Beamsplitter,cuts::Tuple{T,T},idf::DataFrame,cdf::DataFrame,ftdf::DataFrame;update::String = "cat") where {T<:Integer}
    fieldInstance!(efieldp,n,bs,idf,update = update)
    corrInstance!(idf,cdf,update=update)
    ftInstance!(cdf,ftdf,cuts,update=update)
    return idf,cdf,ftdf
end

function run()

    ############################################################################ 
    # Specify parameters
    ############################################################################ 
    makeInstances = true

    tres = 0.01 # 10 picoseconds
    tmax = 10.0 # in nanoseconds
    
    # set simulation seed (-1 for arbitrary seed)
    seed = -1

    # set the number of emitting atoms
    bigN = 50

    # Balmer-α lines specified here
    ωM = [456811.0, 456812.0]
    shift = ωM[2]-ωM[1]
    ωM = 2*π*ωM
    
    # magnitude of each line
    mag = convert(Vector{ComplexF64},ones(length(ωM)))

    # Doppler broadening
    temp = 1000.0 # Kelvins
    σDopp = σTemp(ωM[1],temp)

    # photon counts
    nbar = 10
    ntot = 1000

    # define beamsplitter
    bs = Beamsplitter(0.5,0.5)

    # low cut on τ (indices from start)
    τstart =100
    # high cut on τ (indices from end)
    τend = 0

    ############################################################################ 
    # Set up simulation
    ############################################################################ 

    # Store parameters
    eParams = eFieldParams(mag,ωM,σDopp)
    seed!(eParams,seed)

    # generate times in tres ps intervals up to 2*tmax
    times = collect(0:tres:2*tmax);

    # dataframe for intensity
    dfIntensity = DataFrame(:time=>times)

    # limit the window to tmax to avoid correlation cutoff
    window = convert(Integer,floor(length(times)/2));

    # τ is just the times up to our window
    τ = times[1:window];
    # dataframe for tau dependent values
    dfCorr = DataFrame(:tau=>τ)

    # g2τCalc = map(tau->g2Calc(tau,bigN,eParams),τ)

    # generate frequency bins of fourier transform
    allfreqs = fftFreq(tres,τ,(τstart,τend))
    
    # dataframe for fourier transform results
    dfFreqs = DataFrame(:freq=>allfreqs)

    # concatenation of all photon correlation times
    allCorrTimes = Vector{Float64}[]


    # make first instance of classical calculations   
    classicalInstance!(eParams,bigN,bs,(τstart,τend),dfIntensity,dfCorr,dfFreqs,update = "sum")

	# for n = 1:trials


		# # calculate the average photon counts in each time bin from the beam intensity and the overall average photon count rate
		# γavg1Beam = γIntensity(intensity1Beam,nbar/2)

		# # generate counts for each beam
		# γcounts1Beam1 = poissonCount.(γavg1Beam)
		# γcounts1Beam2 = poissonCount.(γavg1Beam)

		# # concatenate the correlations from this run into all others
		# allCorrTimes = vcat(allCorrTimes,corrTimes(τ,γcounts1Beam1,γcounts1Beam2))

		# # add counts from both beams
		# totCounts .+= γcounts1Beam1
		# totCounts .+= γcounts1Beam2

		# # reinstantiate if desired
		# if makeInstances && n != trials
			# # reinstantiate electric field
			# eInstance1 = eFieldInstance(bigN,eParams)

			# # calculate electric field vs time
			# eFieldT1 = map(t->electricField(t,eInstance1),times)

			# # apply beam splitter
			# eFieldT1Beam = bs.t * eFieldT1

			# # calculate intensity
			# intensity1Beam = intensity.(eFieldT1Beam)

			# # calculate classical g2τ
			# g2τ1Norm = mean(intensity1Beam)^2
			# g2τ1 = map(i->autocorrelate(intensity1Beam,i,window),collect(0:window-1))/g2τ1Norm

			# # calculate Fourier transform of this instance
			# g2τ1FFT = meanFFT(g2τ1,(τstart,τend))

			# # accumulate Fourier transform sum
			# g2τ1FFTsum .+= g2τ1FFT
		# end
	# end
	# # end

	# g2τ1FFTsinglePos,allfreqsPos = fftPositiveFreq(g2τ1FFT,allfreqs)

    # # bin correlation times into a histogram
	# corrHist = fit(Histogram,allCorrTimes,vcat(τ,τ[end]+tres),closed=:left)

	# # normalize the histogram
	# normCorr = corrHist.weights/sum(corrHist.weights)

	# # subtract mean from histogram to reduce constant term
	# normAvgCorr = normCorr .- mean(normCorr)

	# # take the fourier transform
	# γfft = fft(normAvgCorr)

	# # recover fourier transform frequencies
	# γfreqs = fftfreq(length(γfft),1/tres)

	# # select positive frequencies for plotting
	# γfftPos,γfreqsPos = fftPositiveFreq(γfft,γfreqs)


    # γCorrFreqPlot = plot(γfreqsPos,abs.(γfftPos),label = false)
	# xlabel!(γCorrFreqPlot,"frequency (GHz)")
	# title!(γCorrFreqPlot,"Fourier transform of photon correlations")
	# vline!(γCorrFreqPlot,[shift],label="Frequency shift",ls=:dash)


end

export run

end

import .Tst

Tst.run()
