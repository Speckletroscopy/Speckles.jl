### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 92390b64-820e-4866-ac3b-d84392ed8c1f
using Speckles

# ╔═╡ 065faf5b-b1c4-455b-8cf7-22bd5346f32f
νHα2 = [456810,456813] #GHz

# ╔═╡ a0adbbfc-e9da-4793-b6ab-aade13ff0944
paramDict_n = Dict(
				 :n    => [10,20,30,40,50,60,70,80,90,100], # number of atoms
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
				:repeat     => [100], # number of times to repeat correlation measurement
				:reinstance => [true] # control whether or not frequencies and phases should be reinstanced between measurements
				)


# ╔═╡ efe88451-da29-4d70-932c-ab102a98fd04
Speckles.run(paramDict_n)

# ╔═╡ dc7a78bc-ba86-495b-b771-5d1672e3bdc2
paramDict_fγ = Dict(
				 :n    => [50], # number of atoms
				:νm   => [νHα2], # line frequencies in GHz
				:Em   => ["ones"], # relative line magnitudes
				:σ    => [20.0], # Doppler broadening in GHz
				:fγ   => [2e6,2e5,2e4,2e3,2e2,2e1,2.0,0.2],#,"shot10%","shot50%",10.0,1.0,0.16], # mean photon count rate in GHz
				:deadtime   => [0.0], # detector deadtime in nanoseconds
				:resolution => [0.010],#,0.10], # detector resolution in nanoseconds
				:jitter     => [0.015], # detector timing jitter in nanoseconds 
				:efficiency => [0.9], # detector efficiency
				:darkcounts => [1.0e-8], # detector dark count rate in GHz
				:duration   => [20.0], # duration of each correlation measurement in nanoseconds
				:window     => ["halfwindow"], # time over which to average correlations in nanoseconds
				:repeat     => [100], # number of times to repeat correlation measurement
				:reinstance => [true] # control whether or not frequencies and phases should be reinstanced between measurements
				)

# ╔═╡ 301d70aa-137f-48ac-a33d-569c03b3f6cb
Speckles.run(paramDict_fγ)

# ╔═╡ d6fa02d0-ae58-43f5-9deb-a5d99f776ba6
paramDict_res = Dict(
				 :n    => [50], # number of atoms
				:νm   => [νHα2], # line frequencies in GHz
				:Em   => ["ones"], # relative line magnitudes
				:σ    => [20.0], # Doppler broadening in GHz
				:fγ   => [2e6],#,"shot10%","shot50%",10.0,1.0,0.16], # mean photon count rate in GHz
				:deadtime   => [0.0], # detector deadtime in nanoseconds
				:resolution => [0.010,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1],#,0.10], # detector resolution in nanoseconds
				:jitter     => [0.015], # detector timing jitter in nanoseconds 
				:efficiency => [0.9], # detector efficiency
				:darkcounts => [1.0e-8], # detector dark count rate in GHz
				:duration   => [20.0], # duration of each correlation measurement in nanoseconds
				:window     => ["halfwindow"], # time over which to average correlations in nanoseconds
				:repeat     => [100], # number of times to repeat correlation measurement
				:reinstance => [true] # control whether or not frequencies and phases should be reinstanced between measurements
				)

# ╔═╡ b24b54b2-ca4c-442b-8ddf-141aa6a710c8
tbl = Speckles.run(paramDict_res)

# ╔═╡ 50252161-1cca-453a-a12c-5f544b117b1f
begin
	length(SpeckleParamsVector(paramDict_n))+length(SpeckleParamsVector(paramDict_fγ))+length(SpeckleParamsVector(paramDict_res))
end

# ╔═╡ 85f08206-e065-4f76-a9d1-802fb97e0bec
tblNoMissing = dropmissing(tbl)

# ╔═╡ d73899fa-692b-479f-94b3-3efbc79cede2


# ╔═╡ Cell order:
# ╠═92390b64-820e-4866-ac3b-d84392ed8c1f
# ╠═065faf5b-b1c4-455b-8cf7-22bd5346f32f
# ╠═a0adbbfc-e9da-4793-b6ab-aade13ff0944
# ╠═efe88451-da29-4d70-932c-ab102a98fd04
# ╠═dc7a78bc-ba86-495b-b771-5d1672e3bdc2
# ╠═301d70aa-137f-48ac-a33d-569c03b3f6cb
# ╠═d6fa02d0-ae58-43f5-9deb-a5d99f776ba6
# ╠═b24b54b2-ca4c-442b-8ddf-141aa6a710c8
# ╠═50252161-1cca-453a-a12c-5f544b117b1f
# ╠═85f08206-e065-4f76-a9d1-802fb97e0bec
# ╠═d73899fa-692b-479f-94b3-3efbc79cede2
