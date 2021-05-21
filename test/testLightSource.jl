n = 10
EmError = ones(5)
νm = [456812, 456808,456811, 456802]
Em = ones(length(νm))
σ = 20.0
γRate = 1.0e8
source = LightSource(n,Em,νm,σ,γRate)


# test LightSource assertion errors
@test_throws AssertionError LightSource(n,EmError,νm,σ,γRate)
@test_throws AssertionError LightSource(n,Em,νm,σ,456900,456800,γRate)

# test derived quantity functions
@test ν0(source) == 456808.25
@test Δm(source) == [3.75, -0.25, 2.75, -6.25]
@test lineShifts(source) == [4, 1, 10, 3, 6, 9]

# test eField construction
field = eField(source,42)
ν1 = field.νn[1]
ϕ11 = field.ϕmn[1,1]
ν1test = 456797.12946247705
ϕ11test = 5.89027451092873
@test ν1≈ν1test
@test ϕ11≈ϕ11test

@test intensity(0.5) == 0.25
@test intensity(0.01,field) ≈ 6.060748222475236

t = 10.0 #nanoseconds
dt = 0.1 #nanoseconds
nbar = 100.0
intensityVec = map(time->intensity(time,field),0:dt:t)
γVec = γIntensity(nbar,intensityVec)
@test sum(γVec.γvec)≈nbar

bs =  Beamsplitter(1,1)
@test bs.r^2+bs.t^2 ≈ 1.0