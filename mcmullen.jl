### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ ec15aafa-fb76-4004-9102-7cebe68c93d1
using LinearAlgebra

# ╔═╡ 8a20d044-42ac-4614-89aa-a06c4af1e83f
using SparseArrays

# ╔═╡ 4d3a4271-4394-4a94-a516-cff6b3675297
struct Node
	tuple
	children
	word
	valid
end

# ╔═╡ f63f6a13-a315-4887-897a-7559e35d17e9
numbers = Dict()

# ╔═╡ 2139ffed-dba3-433e-ae66-6c15948b36f9
whole_tree = Dict()

# ╔═╡ b8febb46-522e-42da-8d7b-0d54d27249bb
point_list = [0, 0+0.5im, -1*sqrt(3)/4-0.25im, sqrt(3)/4-0.25im]

# ╔═╡ 0c0e9037-99b0-4ea6-adae-110d4ee17b3b
gen1 = stack([[-1.,2.,2.,2.],[0,1.,0,0],[0,0,1.,0],[0,0,0,1.]],dims=1)

# ╔═╡ 23ed120c-452b-4d95-9ff7-ef347bc11c0d
gen2 = stack([[1.,0,0,0],[2.,-1.,2.,2.],[0,0,1.,0],[0,0,0,1.]],dims=1)

# ╔═╡ 1fd4a0f3-2f69-4fcf-87d4-2c940edae987
gen3 = stack([[1.,0,0,0],[0,1.,0,0],[2.,2.,-1.,2.],[0,0,0,1.]],dims=1)

# ╔═╡ 8b8d0b26-5255-46e6-97da-b82faa282128
gen4 = stack([[1.,0,0,0],[0,1.,0,0],[0,0,1.,0],[2.,2.,2.,-1.]],dims=1)

# ╔═╡ 32f127fb-cbf5-4667-92dc-610193e2b655
gen_list = [gen1, gen2, gen3, gen4]

# ╔═╡ 452f4307-9a14-420f-b0c0-cf90f22370b2
function func1(z)
	return conj((2 - sqrt(3))^2 / z)
end

# ╔═╡ 51752925-ecb0-472d-b124-784cb49f6772
function func2(z)
	return 2im + conj(3 / (z - 2im))
end

# ╔═╡ 6c3ef2cd-f08b-448e-a9ee-62c3fb02dd65
function func3(z)
	return -1*sqrt(3) - 1im + conj(3 / (z + sqrt(3) + 1im))
end

# ╔═╡ 0d72c245-2c83-4933-920e-8c5c3ea05dea
function func4(z)
	return sqrt(3) - 1im + conj(3 / (z - sqrt(3) + 1im))
end

# ╔═╡ d0ce22ea-6e8f-4815-84c5-2cb21994bf6c
func_list = [func1, func2, func3, func4]

# ╔═╡ e1f8c861-4545-4694-81f6-a49c757d73dd
function der1(z)
	return abs(z^2 / (2 - sqrt(3))^2)
end

# ╔═╡ 32c1fb88-fc70-46aa-8c11-0e5ba921849c
function der2(z)
	return abs((z - 2im)^2 / 3)
end

# ╔═╡ 9e0b0873-2781-43ee-8f3a-14a69f2d0d03
function der3(z)
	return abs((z + sqrt(3) + 1im)^2 / 3)
end

# ╔═╡ fab98970-cd5e-47cd-a083-ce9c391deb6c
function der4(z)
	return abs((z - sqrt(3) + 1im)^2 / 3)
end

# ╔═╡ c1fb52eb-f5fe-4ec5-b9d8-6b08b3fe7a16
der_list = [der1, der2, der3, der4]

# ╔═╡ c4cb2d62-e407-441c-929f-eef11f3facb6
"maps the starting sample points to their images under the given set of generators"
function sample_point(word)
	p = point_list[last(word)]
    for letter in reverse(word)[2:end]
        p = func_list[letter](p)
	end
    return p
end

# ╔═╡ dcf16831-3f44-4a9c-a9ae-8533630025b4
"determines the value of f'_i(y_{ij})"
function sample_value(word)
	return der_list[word[1]](sample_point(word))
end

# ╔═╡ beffc494-2e0c-4051-9853-729224e3f41b
function secant(x0,y0,x1,y1,z)
	return x0 - (y0-z) * ((x1-x0)/(y1-y0))
end

# ╔═╡ c53c707b-5584-4269-9b59-85dac2bea60b
function power_method(col,row,data,l,a)
	matr = sparse(col, row, data.^a)
    vec = ones(l)
    previous_entry = vec[1]
    previous_val = 0
    current = matr * vec
    current_val = current[1] / previous_entry
    count = 0
    while count < 1000 && abs(current_val - previous_val) > 1e-5
        previous_val = current_val
        previous_entry = current[1]
        current = matr * current
        current_val = current[1] / previous_entry
        count += 1
	end
    return current_val
end

# ╔═╡ 966e4a3e-e84c-4732-9264-22e0302e35a9
function secant_method(col,row,data,l,z,x1,x2,e,its)
	k1 = x1
	k2 = x2
	y1 = power_method(col,row,data,l,k1)
	y2 = power_method(col,row,data,l,k2)
	count = 1
	while abs(y1-z)>e && count<its
		k3 = secant(k1,y1,k2,y2,z)
		k1 = k2
		y1 = y2
		k2 = k3
		y2 = power_method(col,row,data,l,k2)
		count += 1
	end
	return count,k1,y1
end

# ╔═╡ 5080b2f3-5329-40a5-be94-52ae9364a9d2
dual_curvature = 500

# ╔═╡ e214345a-0e18-464c-96f7-30d270903698
"determines if a tuple has too large of a dual curvature"
function kill(tuple,g)
	temp = deleteat!(tuple,g)
	return sqrt(temp[1]*temp[2] + temp[1]*temp[3] + temp[2]*temp[3]) < dual_curvature
end

# ╔═╡ 420ba9aa-2b9a-418c-a7ad-39480bcfb778
"if a node is currently a leaf, attempts to generate its children"
function spawn(leaf)
	if !leaf.valid
		return leaf
	else
		for g in 1:4
			if length(leaf.word) == 0 || g != last(leaf.word)
				child_word = push!(copy(leaf.word),g)
				child = Node(gen_list[g]*leaf.tuple, [], child_word, kill(gen_list[g]*leaf.tuple,g))
				push!(leaf.children, child)
				whole_tree[join(child_word)] = leaf
			end
		end
	end
	return leaf
end

# ╔═╡ fbf7e09c-4916-4bdc-a763-fa5e55b8814e
"generates the tree"
function construct_tree()
	root = Node([-1., 2+sqrt(3), 2+sqrt(3), 2+sqrt(3)],[],[],true)
	current_leaves = [root]
	new_leaves = []
	num_nodes = 1
	while true
		new_leaves = []
		for leaf in current_leaves
			if !leaf.valid
				push!(new_leaves, leaf)
			else
				leaf = spawn(leaf)
				new_leaves = vcat(new_leaves, leaf.children)
				num_nodes += length(leaf.children)
			end
		end
		temp = 0
		for leaf in new_leaves
			if leaf.valid == false
				temp += 1
			end
		end
		if temp == length(new_leaves)
			break
		else
			current_leaves = copy(new_leaves)
		end
	end

	for (i,node) in enumerate(new_leaves)
		numbers[join(node.word)] = i
	end
	return new_leaves
end

# ╔═╡ b5bab224-4e0a-4f0a-9b39-db0737352b8d
"generates the data for the transition matrix"
function construct_data()
	data = construct_tree()
	col = []
	row = []
	nonzero = []
	for (i,node) in enumerate(data)
		sample = join(node.word[2:end])
		for leaf in data
			if startswith(join(leaf.word), sample)
				push!(row, i)
				push!(col, numbers[join(leaf.word)])
				push!(nonzero, sample_value(leaf.word))
			end
		end
	end
	return col, row, nonzero, length(data)
end

# ╔═╡ 0d11aa64-22f5-421f-9ef2-96efa095d689
col, row, data, l = construct_data()

# ╔═╡ bf8ba11c-c87a-42e1-8b79-590b513959a7
secant_method(col,row,data,l,1,1.3,1.31,1e-7,10000)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.4"
manifest_format = "2.0"
project_hash = "6865b0a5e609d4d4358ebe0134f60aa2e67c548e"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═ec15aafa-fb76-4004-9102-7cebe68c93d1
# ╠═8a20d044-42ac-4614-89aa-a06c4af1e83f
# ╠═4d3a4271-4394-4a94-a516-cff6b3675297
# ╠═f63f6a13-a315-4887-897a-7559e35d17e9
# ╠═2139ffed-dba3-433e-ae66-6c15948b36f9
# ╠═d0ce22ea-6e8f-4815-84c5-2cb21994bf6c
# ╠═c1fb52eb-f5fe-4ec5-b9d8-6b08b3fe7a16
# ╠═32f127fb-cbf5-4667-92dc-610193e2b655
# ╠═b8febb46-522e-42da-8d7b-0d54d27249bb
# ╠═0c0e9037-99b0-4ea6-adae-110d4ee17b3b
# ╠═23ed120c-452b-4d95-9ff7-ef347bc11c0d
# ╠═1fd4a0f3-2f69-4fcf-87d4-2c940edae987
# ╠═8b8d0b26-5255-46e6-97da-b82faa282128
# ╠═452f4307-9a14-420f-b0c0-cf90f22370b2
# ╠═51752925-ecb0-472d-b124-784cb49f6772
# ╠═6c3ef2cd-f08b-448e-a9ee-62c3fb02dd65
# ╠═0d72c245-2c83-4933-920e-8c5c3ea05dea
# ╠═e1f8c861-4545-4694-81f6-a49c757d73dd
# ╠═32c1fb88-fc70-46aa-8c11-0e5ba921849c
# ╠═9e0b0873-2781-43ee-8f3a-14a69f2d0d03
# ╠═fab98970-cd5e-47cd-a083-ce9c391deb6c
# ╠═c4cb2d62-e407-441c-929f-eef11f3facb6
# ╠═dcf16831-3f44-4a9c-a9ae-8533630025b4
# ╠═e214345a-0e18-464c-96f7-30d270903698
# ╠═420ba9aa-2b9a-418c-a7ad-39480bcfb778
# ╠═fbf7e09c-4916-4bdc-a763-fa5e55b8814e
# ╠═b5bab224-4e0a-4f0a-9b39-db0737352b8d
# ╠═beffc494-2e0c-4051-9853-729224e3f41b
# ╠═c53c707b-5584-4269-9b59-85dac2bea60b
# ╠═966e4a3e-e84c-4732-9264-22e0302e35a9
# ╠═5080b2f3-5329-40a5-be94-52ae9364a9d2
# ╠═0d11aa64-22f5-421f-9ef2-96efa095d689
# ╠═bf8ba11c-c87a-42e1-8b79-590b513959a7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
