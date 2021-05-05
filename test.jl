module Tst

using Speckles
using Plots

function run()

    ############################################################################ 
    # Specify parameters
    ############################################################################ 
    makeInstances = true

    tres = 0.01 # 10 picoseconds
    tmax = 20.0 # in nanoseconds
    
    # set simulation seed (-1 for arbitrary seed)
    seed = -1

    # set the number of emitting atoms
    bigN = 100

    # Balmer-α lines specified here
    ωM = [456811.0, 456815.0]
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

    # have enough trials to get the desired number of photon counts
    trials = convert(Integer,ceil(ntot/nbar))

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

    # array to accumulate all photon counts
	totCounts = zeros(Integer,length(times))

    # concatenation of all photon correlation times
    allCorrTimes = Vector{Float64}[]

    ud = "sum"
    # make first instance of classical calculations   
    classicalInstance!(eParams,bigN,bs,(τstart,τend),dfIntensity,dfCorr,dfFreqs,update = ud)

    for n = 1:trials

        # calculate the average photon counts in each time bin from the beam intensity and the overall average photon count rate
        γavgBeam = γIntensity(dfIntensity[!,end],nbar/2)

        # generate counts for each beam
        γcountsBeam1 = poissonCount.(γavgBeam)
        γcountsBeam2 = poissonCount.(γavgBeam)

        # concatenate the correlations from this run into all others
        allCorrTimes = vcat(allCorrTimes,corrTimes(τ,γcountsBeam1,γcountsBeam2))

        # add counts from both beams
        totCounts .+= γcountsBeam1
        totCounts .+= γcountsBeam2

        # reinstantiate if desired
        if makeInstances && n != trials
            classicalInstance!(eParams,bigN,bs,(τstart,τend),dfIntensity,dfCorr,dfFreqs,update = ud)
        end
    end
    

    ftSum, ftFreqs = fftPositiveFreq(dfFreqs.sum,dfFreqs.freq)
    sumplot = plot(ftFreqs,abs.(ftSum))
    savefig(sumplot,"sumplot.svg")

    single = dfFreqs.PowerSpec
    ftSingle, ftFreqs = fftPositiveFreq(single,dfFreqs.freq)
    singleplot = plot(ftFreqs,abs.(ftSingle))
    savefig(singleplot,"singleplot.svg")

end

export run

end

import .Tst

Tst.run()
