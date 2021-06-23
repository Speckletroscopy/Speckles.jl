### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 736da374-bfee-11eb-23a7-f1ad290a242e
using Speckles

# ╔═╡ feeaa93e-1a23-4cfc-b219-701531191d10
begin
	# specify all parameters
	νHα1 = [456808] #GHz
	νHα2 = [456808,456811] #GHz
	νHα = [456812, 456808,456811, 456802] #GHz

	paramDict = Dict(
					:n    => [10], # number of atoms
					:νm   => [νHα1,νHα2], # line frequencies in GHz
					:Em   => ["ones"], # relative line magnitudes
					:σ    => [20.0], # Doppler broadening in GHz
					:fγ   => [2*1e4/0.01], # mean photon count rate in GHz
					:deadtime   => [0.00],#,0.01,0.025,0.05,0.1], # detector deadtime in nanoseconds
					:resolution => [0.010], # detector resolution in nanoseconds
					:jitter     => [0.015], # detector timing jitter in nanoseconds 
					:efficiency => [0.9], # detector efficiency
					:darkcounts => [1.0e-8], # detector dark count rate in GHz
					:duration   => [20.0], # duration of each correlation measurement in nanoseconds
					:repeat     => [10], # number of times to repeat correlation measurement
					:reinstance => [true], # control whether or not frequencies and phases should be reinstanced between measurements
					:timeint    => ["halfwindow"] # time over which to average correlations in nanoseconds
					# :directory  => [] # defaults to main package directory
					)

end

# ╔═╡ 5e100134-d913-44a5-ad2b-55d29e889304
# split into vector of dictionaries: one for each run
iterParams = paramVector(paramDict);

# ╔═╡ 4d9ab06d-6bae-4066-9333-4287a948112e
# define beamsplitter
bs = Beamsplitter(0.5,0.5)

# ╔═╡ adc32c75-b46f-4818-8de8-50d9990b26b0
for i=1:length(iterParams)
	iterParams[i][:Em] = ones(length(iterParams[i][:νm]))
	iterParams[i][:timeint] = iterParams[i][:duration]/2
end

# ╔═╡ d237c492-4bc3-45a7-a3b5-63a33988a94f
prefixes = map(x->makeName(x),iterParams)

# ╔═╡ e893f820-4be5-4810-84a9-b5629119fe1f
function makeSource(params::Dict)
	return LightSource(
				params[:n],
				params[:Em],
				params[:νm],
				params[:σ],
				params[:fγ]
			)
end

# ╔═╡ 4a19ac54-9617-4a29-a9e6-04cb91450470
source = makeSource(iterParams[1])

# ╔═╡ f4f1c6d6-2b3a-41b0-88e5-d0a59debed04
function makeDetector(params::Dict)
	return Detector(
            params[:deadtime],
            params[:resolution],
            params[:jitter],
            params[:efficiency],
            params[:darkcounts]
        )
end

# ╔═╡ dac6bbdc-8996-425b-8410-f789b71ea3d9
detect = makeDetector(iterParams[1])

# ╔═╡ 65bf052e-43c6-4743-b264-b04ae9843a2c
γint = γIntensity(nbar(iterParams[1][:duration],source)*bs.t,
	iterParams[1][:duration],
	detect.resolution,
	source)

# ╔═╡ b10ddbe6-90eb-4e65-a972-d20f52ed7547


# ╔═╡ Cell order:
# ╠═736da374-bfee-11eb-23a7-f1ad290a242e
# ╠═feeaa93e-1a23-4cfc-b219-701531191d10
# ╠═5e100134-d913-44a5-ad2b-55d29e889304
# ╠═4d9ab06d-6bae-4066-9333-4287a948112e
# ╠═adc32c75-b46f-4818-8de8-50d9990b26b0
# ╠═d237c492-4bc3-45a7-a3b5-63a33988a94f
# ╠═e893f820-4be5-4810-84a9-b5629119fe1f
# ╠═4a19ac54-9617-4a29-a9e6-04cb91450470
# ╠═f4f1c6d6-2b3a-41b0-88e5-d0a59debed04
# ╠═dac6bbdc-8996-425b-8410-f789b71ea3d9
# ╠═65bf052e-43c6-4743-b264-b04ae9843a2c
# ╠═b10ddbe6-90eb-4e65-a972-d20f52ed7547
