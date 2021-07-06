νHα2 = [456810,456813] #GHz
paramDict = Dict(
                :n    => [100], # number of atoms
                :νm   => [νHα2], # line frequencies in GHz
                :Em   => ["ones"], # relative line magnitudes
                :σ    => [20.0], # Doppler broadening in GHz
                :fγ   => [2.0e6],#,"shot10%","shot50%",10.0,1.0,0.16], # mean photon count rate in GHz
                :deadtime   => [0.0], # detector deadtime in nanoseconds
                :resolution => [0.010],#,0.10], # detector resolution in nanoseconds
                :jitter     => [0.015], # detector timing jitter in nanoseconds 
                :efficiency => [0.9], # detector efficiency
                :darkcounts => [1.0e-8], # detector dark count rate in GHz
                :duration   => [20.0], # duration of each correlation measurement in nanoseconds
                :window     => ["halfwindow"], # time over which to average correlations in nanoseconds
                :repeat     => collect(10:10:100), # number of times to repeat correlation measurement
                :reinstance => [true] # control whether or not frequencies and phases should be reinstanced between measurements
                )
nameDict = Dict(
                :n    => "number of atoms",
                :νm   => "line frequencies (GHz)",
                :fγ   => "mean photon count rate (GHz)",
                :deadtime   => "dead time (ns)",
                :resolution => "resolution (ns)",
                :duration   => "measurement duration (ns)",
                :repeat     => "trials",
                :reinstance => "randomize each trial"
               )
