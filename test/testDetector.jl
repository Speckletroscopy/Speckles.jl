
# light source parameters
n = 10
νm = [456812, 456808,456811, 456802]
Em = ones(length(νm))
σ = 20.0
γRate = 1.0e8


# detector parameters
deadtime = 10.0 #nanoseconds
resolution = 0.01 #nanoseconds
jitter = 0.015 #nanoseconds
efficiency = 0.9
darkcounts = 1.0e-8 #GHz

# Detector object construction
@test_throws AssertionError Detector(deadtime,resolution,jitter,1.1,darkcounts)

# create detector and light source objects
source = LightSource(n,Em,νm,σ,γRate)
detect = Detector(deadtime,resolution,jitter,efficiency,darkcounts)

readout(20,source,detect)