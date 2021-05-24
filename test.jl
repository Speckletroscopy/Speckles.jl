module Tst

using Speckles

function sigmaTemp5k()
    σTemp(456811.0,5000.0)
end
export sigmaTemp5k

function run()
    # specify all parameters
    νHα1 = [456808] #GHz
    νHα2 = [456808,456811] #GHz
    # νHα = [456812, 456808,456811, 456802] #GHz

    paramDict = Dict(
                    :n    => [10], # number of atoms
                    :νm   => [νHα2], # line frequencies in GHz
                    :Em   => ["ones"], # relative line magnitudes
                    :σ    => [20.0], # Doppler broadening in GHz
                    :fγ   => ["shot1%"], # mean photon count rate in GHz
                    :deadtime   => [0.00,0.01,0.025,0.05,0.1], # detector deadtime in nanoseconds
                    :resolution => [0.010], # detector resolution in nanoseconds
                    :jitter     => [0.015], # detector timing jitter in nanoseconds 
                    :efficiency => [0.9], # detector efficiency
                    :darkcounts => [1.0e-8], # detector dark count rate in GHz
                    :duration   => [20.0], # duration of each correlation measurement in nanoseconds
                    :repeat     => [10], # number of times to repeat correlation measurement
                    :reinstance => [true], # control whether or not frequencies and phases should be reinstanced between measurements
                    :timeint    => ["halfwindow"] # time over which to average correlations in nanoseconds
                    )

    # split into vector of dictionaries: one for each run
    iterParams = paramVector(paramDict)
    
    # define beamsplitter
    bs = Beamsplitter(0.5,0.5)

    # iterate over parameter dictionaries
    for params in iterParams
        if params[:Em] == "ones" # sets all line magnitudes equal to one if true
            params[:Em] = ones(length(params[:νm]))
        end
        if params[:fγ] == "shot1%"
            params[:fγ] = 2*1e4/params[:resolution] # multiply by 2 so error in each beam is ~1%
        elseif params[:fγ] == "shot10%"
            params[:fγ] = 2*1e2/params[:resolution] # multiply by 2 so error in each beam is ~10%
        elseif params[:fγ] == "shot50%"
            params[:fγ] = 2*4/params[:resolution] # multiply by 2 so error in each beam is ~50%
        end
        if params[:timeint] == "halfwindow"
           params[:timeint] = params[:duration]/2 
        end

        # create naming prefix
        prefix = makeName(params)

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

        # calculate the length in indices of the time integration window
        indexint = length(γint)÷convert(Int,ceil(params[:duration]/params[:timeint]))

        # generate the correlation time offsets
        τ = params[:resolution]*collect(0:(indexint-1))
        
        # create a 5 e-fold low cut in τ to get rid of non-periodic part of g2
        τCutLow = sqrt(5)/params[:σ] # low τ cut in nanoseconds

        # calculate index of τ cut
        iτCutLow = convert(Int,ceil(τCutLow/params[:resolution]))

        # get the frequencies for a Fourier transform over τ with cuts
        freqs = fftFreq(params[:resolution],τ,(iτCutLow,lastindex(τ)))

        # initialize time and frequency DataFrames
        timeData = DataFrame(:time=>τ)
        freqData = DataFrame(:freq=>freqs[freqs .>= 0])
        
        # iterate over desired number of repeats
        for i=1:params[:repeat]
            # read out each beam
            readout1 = denseReadout(γint,detect)
            readout2 = denseReadout(γint,detect)
            # calculate correlation
            corr12 = map(offset->ncorrelate(readout1,readout2,offset,indexint),0:(indexint-1))
            # store instance in dataframe
            iname = "corr$i"
            timeData[!,iname] = corr12
            
            # take the Fourier transform of the correlation function
            ftcorr12 = meanFFT(corr12,(iτCutLow,lastindex(corr12)))
            ftcorr12pos, = fftPositiveFreq(ftcorr12,freqs) 
            
            freqData[!,iname] = abs.(ftcorr12pos)

            # generate new intensity instance if desired
            if params[:reinstance]
                γint = γIntensity(
                    nbar(params[:duration],source)*0.5, # multiply by 0.5 for two beams 
                    params[:duration],
                    detect.resolution,
                    source
                    )
            end
        end
        
        # create sum column in frequency and time DataFrames
        timeData[!,:sum] = sum(eachcol(timeData)[2:end])
        freqData[!,:sum] = sum(eachcol(freqData)[2:end])
        
        # make names for data files
        timeDataName = dataDir(string(prefix,"g2-vs-tau.csv"))
        freqDataName = dataDir(string(prefix,"ftg2-vs-freq.csv"))
        
        # save data
        CSV.write(timeDataName,timeData)
        CSV.write(freqDataName,freqData)
        
        # create and save plots
        γCorrTimePlot(timeData, params, prefix)
        γCorrFreqPlot(freqData, params, prefix)
    end
    return nothing
end #run
export run
end #Tst
import .Tst

Tst.run()