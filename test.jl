module Tst

using Speckles

function sigmaTemp5k()
    σTemp(456811.0,5000.0)
end
export sigmaTemp5k

function run3()
    # specify all parameters
    νHα = [456808,456811] #GHz
    # νHα = [456812, 456808,456811, 456802] #GHz
    EmHα = ones(length(νHα))

    paramDict = Dict(
                    :n    => [10], # number of atoms
                    :νm   => [νHα], # line frequencies in GHz
                    :Em   => [EmHα], # relative line magnitudes
                    :σ    => [20.0], # Doppler broadening in GHz
                    :fγ   => ["shot1%"], # mean photon count rate in GHz
                    :deadtime   => [0.00,0.01,0.05], # detector deadtime in nanoseconds
                    :resolution => [0.010], # detector resolution in nanoseconds
                    :jitter     => [0.015], # detector timing jitter in nanoseconds 
                    :efficiency => [0.9], # detector efficiency
                    :darkcounts => [1.0e-8], # detector dark count rate in GHz
                    :duration   => [20.0], # duration of each correlation measurement in nanoseconds
                    :repeat     => [1], # number of times to repeat correlation measurement
                    :reinstance => [false], # control whether or not frequencies and phases should be reinstanced between measurements
                    :timeint    => ["halfwindow"] # time over which to average correlations in nanoseconds
                    )

    # split into vector of dictionaries: one for each run
    iterParams = paramVector(paramDict)
    
    # define beamsplitter
    bs = Beamsplitter(0.5,0.5)

    out = []
    # iterate over parameter dictionaries
    for params in iterParams
        if params[:fγ] == "shot1%"
            params[:fγ] = 2*1e4/params[:resolution] # multiply by 2 so error in each beam is ~1%
        elseif params[:fγ] == "shot10%"
            params[:fγ] = 2*1e2/params[:resolution] # multiply by 2 so error in each beam is ~1%
        elseif params[:fγ] == "shot50%"
            params[:fγ] = 2*4/params[:resolution] # multiply by 2 so error in each beam is ~1%
        end
        if params[:timeint] == "halfwindow"
           params[:timeint] = params[:duration]/2 
        end

        # create LightSource and Detector objects
        source = LightSource(
            params[:n],
            params[:Em],
            params[:νm],
            params[:σ],
            params[:fγ]
        )
        detect = Detector(
            params[:deadtime],
            params[:resolution],
            params[:jitter],
            params[:efficiency],
            params[:darkcounts]
        )

        # calculate the average photon counts for each time bin
        γint = γIntensity(
            nbar(params[:duration],source)*0.5, # multiply by 0.5 for two beams 
            params[:duration],
            detect.resolution,
            source
            )
        
        indexint = length(γint)÷convert(Int,ceil(params[:duration]/params[:timeint]))
        times = params[:resolution]*collect(0:(indexint-1))
        freqs = fftFreq(params[:resolution],times,(1,lastindex(times)))
        timeData = DataFrame(:time=>times)
        freqData = DataFrame(:freq=>freqs)
        for i=1:params[:repeat]
            # read out each beam
            readout1 = denseReadout(γint,detect)
            readout2 = denseReadout(γint,detect)
            corr12 = map(offset->ncorrelate(readout1,readout2,offset,indexint),0:indexint)
            iname = "corr$i"
            timeData[!,iname] = corr12
        end

    end
    return timedata
end

export run3

function run2()

    # specify all parameters
    νHα = [456812, 456808,456811, 456802] #GHz
    EmHα = ones(length(νHα))

    paramDict = Dict(
                    :n    => [10,100], # number of atoms
                    :νm   => [νHα], # line frequencies in GHz
                    :Em   => [EmHα], # relative line magnitudes
                    :σ    => [20.0], # Doppler broadening in GHz
                    :fγ   => [0.16], # mean photon count rate in GHz
                    :deadtime   => [10.0], # detector deadtime in nanoseconds
                    :resolution => [0.010,0.1], # detector resolution in nanoseconds
                    :jitter     => [0.015], # detector timing jitter in nanoseconds 
                    :efficiency => [0.9], # detector efficiency
                    :darkcounts => [1.0e-8], # detector dark count rate in GHz
                    :duration   => [1.0e3], # duration of each correlation measurement in nanoseconds
                    :repeat     => [1], # number of times to repeat correlation measurement
                    :reinstance => [false], # control whether or not frequencies and phases should be reinstanced between measurements
                    :timeint    => [5.0e2] # time over which to average correlations in nanoseconds
                    )

    # split into vector of dictionaries: one for each run
    iterParams = paramVector(paramDict)
    
    # define beamsplitter
    bs = Beamsplitter(0.5,0.5)

    params = iterParams[1]
    # loop over run dictionaries
    # for params in iterParams
        # create LightSource and Detector objects
        source = LightSource(
            params[:n],
            params[:Em],
            params[:νm],
            params[:σ],
            params[:fγ]
        )
        detect = Detector(
            params[:deadtime],
            params[:resolution],
            params[:jitter],
            params[:efficiency],
            params[:darkcounts]
        )

        γint = γIntensity(
            nbar(params[:duration],source)*0.5, # multiply by 0.5 for two beams 
            params[:duration],
            detect.resolution,
            source
            )

        # for i=1:params[:repeat]
        #     if i % 10 == 0
        #         println(i)
        #     end
            # get photon readouts for each detector
            readout1 = readout(γint,detect)
            readout2 = readout(γint,detect)

            return readout1,readout2
            # nzcorrOffset1 = readout2.nzind[1] - readout1.nzind[1]
            # readout1[readout1.nzind[1]]

            # function nzrange(istart::int,irange::int,vec::sparsevector)
            #     nzstart = vec.nzind[istart]
            #     nzend = nzstart+irange
            #     if nzend > length(vec)
            #         nzend = length(vec)
            #     end
            #     return nzstart:nzend
            # end
            
            # function correlatesparsevectors(vec1::sparsevector, vec2::sparsevector)
            #     return sum(vec1 .* vec2)/(sum(vec1)*sum(vec2))
            # end
            # indexint = length(γint)÷convert(int,ceil(params[:duration]/params[:timeint]))
            
            # countvec1 = readout1[nzrange(1,indexint,readout1)]
            # countvec2 = readout2[nzrange(1,indexint,readout2)]


            # correlateSparseVectors(countVec1,countVec2)
            # reinstantiate avg photon intensity if desired
            # if params[:reinstance]
            #     γint = γIntensity(
            #         nbar(params[:duration],source)*0.5,
            #         params[:duration],
            #         detect.resolution,
            #         source
            #         )
            # end
        # end
    # end
    

end

function run()

    ############################################################################ 
    # Specify parameters
    ############################################################################ 
    makeInstances = true

    paramDict = Dict(
                     :tres=>[0.01], # in nanoseconds
                     :tmax=>[10.0], # nanoseconds
                     :bigN=>[10], # number of atoms
                     :mag=>[convert(Vector{ComplexF64},ones(2))], # field magnitude
                     :νM=>[[456811.0,456815.0]], # line frequencies
                     :temp=>[5000], # Kelvins
                     :nbar=>[50], # average photon counts in measurement period
                     :ntot=>[1000], # approximate total photon counts
                     :reset=>[1.0], # detector reset time in nanoseconds
                     :seed=>[-1] # set seed for reproducible results
                    )
    # split paramDict into iterable vector of dictionaries
    iterParams = paramVector(paramDict)

    # define beamsplitter
    bs = Beamsplitter(0.5,0.5)

    # low cut on τ (indices from start)
    τstart =100
    # high cut on τ (indices from end)
    τend = 0

    # set array storage mode
    ud = "sum"

    println("Beginning $(length(iterParams)) simulation runs")
    for (i,params) in enumerate(iterParams)
        ############################################################################ 
        # Set up simulation
        ############################################################################ 
        println("----- Starting run $i -----")

        ωM = 2*π*params[:νM]
        # Calculate thermal Doppler broadening
        σDopp = σTemp(ωM[1],params[:temp])
        # Store parameters
        eParams = eFieldParams(params[:mag],ωM,σDopp)
        seed!(eParams,params[:seed])

        # generate times in tres ps intervals up to 2*tmax
        times = collect(0:params[:tres]:2*params[:tmax]);

        # have enough trials to get the desired number of photon counts
        instances = convert(Integer,ceil(params[:ntot]/params[:nbar]))

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
        allfreqs = fftFreq(params[:tres],τ,(τstart,τend))
        
        # dataframe for fourier transform results
        dfFreqs = DataFrame(:freq=>allfreqs)

        # array to accumulate all photon counts
        totCounts = zeros(Integer,length(times))

        # concatenation of all photon correlation times
        allCorrTimes = Float64[]

        # coincident photon count times
        coincidentTimes = Float64[]

        # make first instance of classical calculations   
        classicalInstance!(eParams,
                           params[:bigN],
                           bs,
                           (τstart,τend),
                           dfIntensity,
                           dfCorr,
                           dfFreqs,
                           update = ud
                          )

        maxTau = 1
        for n = 1:instances

            # calculate the average photon counts in each time bin from the beam intensity and the overall average photon count rate
            γavgBeam = γIntensity(dfIntensity[!,end],params[:nbar]/2)

            # generate counts for each beam
            γcountsBeam1 = poissonCount.(γavgBeam)
            γcountsBeam2 = poissonCount.(γavgBeam)

            γ1 = firstNonzero(γcountsBeam1)
            γ2 = firstNonzero(γcountsBeam2)
            γi = abs(γ2-γ1)+1
            maxTau = γi > maxTau ? γi : maxTau
            # println("Coincident time = ",times[γi])
            push!(coincidentTimes,times[γi])

            # concatenate the correlations from this run into all others
            allCorrTimes = vcat(allCorrTimes,corrTimes(τ,γcountsBeam1,γcountsBeam2))

            # add counts from both beams
            totCounts .+= γcountsBeam1
            totCounts .+= γcountsBeam2

            # reinstantiate if desired
            if makeInstances && n != instances
                classicalInstance!(eParams,params[:bigN],bs,(τstart,τend),dfIntensity,dfCorr,dfFreqs,update = ud)
            end
        end
        
        ########################################################################
        # photon counting
        ########################################################################
        corrHist = fit(Histogram,allCorrTimes,vcat(τ,τ[end]+params[:tres]),closed=:left)
        coincidentHist = fit(Histogram,coincidentTimes,times[1:maxTau+1],closed=:left)

        # normalize the histogram
        normCorr = corrHist.weights/sum(corrHist.weights)

        # subtract mean from histogram to reduce constant term
        normAvgCorr = normCorr .- mean(normCorr)

        # take the fourier transform
        γfft = fft(normAvgCorr)

        # recover fourier transform frequencies
        γfreqs = fftfreq(length(γfft),1/params[:tres])

        # select positive frequencies for plotting
        γfftPos,γfreqsPos = fftPositiveFreq(γfft,γfreqs)

        ########################################################################
        # ---- Make plots and export data ----
        ########################################################################
        prefix = makeName(params)

        ########################################################################
        # classical results
        ########################################################################

        # classical data
        classicalFreqData(dfFreqs,params,prefix)

        # single fourier transform plot
        classicalSinglePlot(dfFreqs,params,prefix)

        # sum of all fourier transforms plot
        classicalSumPlot(dfFreqs,params,prefix)

        ########################################################################
        # photon counting results
        ########################################################################
        γCorrFreqPlot(γfreqsPos,γfftPos,params,prefix)

        coincidentPlot = plot(coincidentHist, label=false)
        xlabel!(coincidentPlot,"Time between counts (ns)")
        title!(coincidentPlot,"Coincident Photon Count Timing")
        plotname = string(prefix,"coincident-counts-vs-time.svg")
        savefig(coincidentPlot,plotsDir(plotname))
        println("----- Finished run $i -----")
    end

end

export run

end

import .Tst
# Tst.run()
# Tst.sigmaTemp5k()
# readout1,readout2 = Tst.run2()
testvec = Tst.run3()
