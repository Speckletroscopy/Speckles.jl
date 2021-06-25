### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 0d1ca7ea-24e8-477d-b0a3-f039641bf24c
include("src/Speckles.jl")

# ╔═╡ 736da374-bfee-11eb-23a7-f1ad290a242e
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end


# ╔═╡ b10ddbe6-90eb-4e65-a972-d20f52ed7547
speckles = ingredients("src/Speckles.jl")

# ╔═╡ 92390b64-820e-4866-ac3b-d84392ed8c1f
speckles.run()

# ╔═╡ Cell order:
# ╠═0d1ca7ea-24e8-477d-b0a3-f039641bf24c
# ╠═736da374-bfee-11eb-23a7-f1ad290a242e
# ╠═b10ddbe6-90eb-4e65-a972-d20f52ed7547
# ╠═92390b64-820e-4866-ac3b-d84392ed8c1f
