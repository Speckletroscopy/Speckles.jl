module Tst

using Speckles

function run()

    ############################################################################ 
    # Specify parameters
    ############################################################################ 
    makeInstances = true
    # paramDict = Dict(
                     # :tres=>[0.01], # in nanoseconds
                     # :tmax=>[10.0,100.0,1000.0], # nanoseconds
                     # :bigN=>[10,50,100,500], # number of atoms
                     # :mag=>[convert(Vector{ComplexF64},ones(3))], # field magnitude
                     # :νM=>[[456811.0,456813.0,456815.0]], # line frequencies
                     # :temp=>[5000], # Kelvins
                     # :nbar=>[10,50,100,500], # average photon counts in measurement period
                     # :ntot=>[10000], # approximate total photon counts
                     # :reset=>[1.0], # detector reset time in nanoseconds
                     # :seed=>[-1] # set seed for reproducible results
                    # )

    paramDict = Dict(
                     :tres=>[0.01], # in nanoseconds
                     :tmax=>[10.0], # nanoseconds
                     :bigN=>[10], # number of atoms
                     :mag=>[convert(Vector{ComplexF64},ones(2))], # field magnitude
                     :νM=>[[456811.0,456815.0]], # line frequencies
                     :temp=>[5000], # Kelvins
                     :nbar=>[10,50,100,500], # average photon counts in measurement period
                     :ntot=>[100000], # approximate total photon counts
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

Tst.run()
