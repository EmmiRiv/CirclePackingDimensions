### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 236bdd31-4496-48e0-be58-7319738fc1b6
using LinearAlgebra

# ╔═╡ 6be3c15b-e874-48b2-9d98-e7b8c2f2621a
using SpecialFunctions

# ╔═╡ 9144b032-ed27-4910-bf0e-b26f1bf3b8ff
r_mat = (1/12)*[[-6+8im,11im],[4im,-6-8im]]

# ╔═╡ 6a4be131-7112-4a03-b30a-e6d1fe43e0f1
r_mat2 = [[-(7 + 24im)/36, -121/144], [-1/9, (-7 + 24im)/36]]

# ╔═╡ 1041df01-3192-4f28-a0d9-77733cbb6139
g_mat = (1/12)*[[2-2im,-1-5im],[4-4im,-2-10im]]

# ╔═╡ aa71d4ec-94aa-4b43-a6cb-1361948d9dc7
Nc = 9

# ╔═╡ 6a94dbda-1ee4-484c-a7cf-2c58aba99ef5
Lc = 20

# ╔═╡ 84fcb6c5-dcb7-49d4-95d8-33aabe5a5a55
k0 = 100

# ╔═╡ 5afdaf31-2bfe-4b1f-8991-a36e1dba5811
print(secant())

# ╔═╡ c8f312cf-b086-4307-8a06-7ff3ef762a56
"Computes the binomial coefficient with a real-valued r"
function f_binomial(r::Real,k::Integer)
	if k<0
		return 0
	end
	res = 1
	for i in 0:(k-1)
		res = res*(r-i)
	end
	return res/factorial(k)
end

# ╔═╡ 77576ece-4e29-46e5-8dc2-35d87ee74449
"Computes a coefficent for the sth term of the power series expansion"
function coefficientM(n::Integer,matr,q::Real,s::Integer)
	total = 0
	for i in 0:n
		addn = binomial(n,i)*f_binomial(-n-q,s-i)*(^(matr[1][1],i))*(^(matr[1][2],n-i))*(^(matr[2][1],s-i))*(^(matr[2][2],-n-q-s+i))
		total = total + addn
	end
	return total
end

# ╔═╡ 0c77f33d-3d6b-4fba-b80e-a5a85460c973
"Creates the list of power series coefficients for the image of one input z*^mz^n"
function generate_func_row(matr,q::Real,n::Integer,m::Integer)
    lst = Vector{ComplexF64}()
    for r_counter in 0:Nc
		for s_counter in 0:Nc
			push!(lst,conj(coefficientM(m,matr,q,r_counter))*coefficientM(n,matr,q,s_counter))
		end
	end
    return lst
end

# ╔═╡ 626dbef0-bc4c-49f1-9eb0-fe9a08fac272
"Generates L_{A1^kR,q}(m,n;r,s) for one k"
function gen_matrix_AkR(matr,q::Real)
	L_k = Vector{Vector{ComplexF64}}()
	for m_counter in 0:Nc
		for n_counter in 0:Nc
			push!(L_k,generate_func_row(matr,q,n_counter,m_counter))
		end
	end
	return L_k
end

# ╔═╡ 746ea648-f98d-47b9-99e8-9655d55186f7
"Generates the sum of the L matrices for k up to k_0"
function gen_approx_k0(q::Real,p::Integer)
	k = 1
	if p == 1
		L = gen_matrix_AkR(r_mat+k*g_mat,q)
		for k in 2:(k0-1)
			L = L + gen_matrix_AkR(r_mat+k*g_mat,q)
		end
	else
		L = gen_matrix_AkR(r_mat2+k*g_mat,q)
		for k in 2:(k0-1)
			L = L + gen_matrix_AkR(r_mat2+k*g_mat,q)
		end
	end
	return L
end

# ╔═╡ c03fd4c5-ad48-4194-9015-ef44b556985f
"Generates the power series coefficent for z^s/k^{q+l}"
function coefficientF(q::Real,l_input::Integer,s::Integer,n::Integer,r_use)
	total = 0
	for i in 0:s
		total_s = 0
		for a in 0:i
			for b in 0:(n-i)
				for c in 0:(s-i)
					d = l_input-a-b-c
					total_s = total_s + binomial(i,a)*binomial(n-i,b)*binomial(s-i,c)*f_binomial(-n-q-s+i,d)*(^(r_use[1][1],a))*(^(g_mat[1][1],i-a))*(^(r_use[1][2],b))*(^(g_mat[1][2],n-i-b))*(^(r_use[2][1],c))*(^(g_mat[2][1],s-i-c))*(^(r_use[2][2],d))*(^(g_mat[2][2],-n-q-s+i-d))
				end
			end
		end
		total = total+binomial(n,i)*f_binomial(-n-q,s-i)*total_s
	end	
	return total
end

# ╔═╡ f1bcffef-817d-44e2-ab36-fc95dfe867e1
"Generates the image of one z*^mz^n for one l"
function generate_func_row2(q::Real,n::Integer,m::Integer,l::Integer,r_use)
	lst = Vector{ComplexF64}()
    for r_counter in 0:Nc
		for s_counter in 0:Nc
			temp_l = 0
			for lp in 0:l
				temp_l = temp_l + conj(coefficientF(q,l-lp,r_counter,m,r_use))*coefficientF(q,l,s_counter,n,r_use)
			end
			push!(lst,temp_l)
		end
	end
    return zeta(l+2*q,k0)*lst
end

# ╔═╡ e66b214d-d997-4afd-b140-062b7938ad7e
"Generates the matrix of power series coefficents for one fixed l"
function gen_matrix_l(q::Real,l::Integer,r_use)
	F_l = Vector{Vector{ComplexF64}}()
	for m_counter in 0:Nc
		for n_counter in 0:Nc
			push!(F_l,generate_func_row2(q,n_counter,m_counter,l,r_use))
		end
	end
	return F_l
end

# ╔═╡ 609b7838-dd9d-4668-af75-25c85d3c0ced
"Generates the second matrix of coefficients"
function gen_approx_F(q::Real,p::Integer)
	l = 0
	if p == 1
		F = gen_matrix_l(q,0,r_mat)
		for l in 1:Lc
			F = F + gen_matrix_l(q,l,r_mat)
		end
	else
		F = gen_matrix_l(q,0,r_mat2)
		for l in 1:Lc
			F = F + gen_matrix_l(q,l,r_mat2)
		end
	end
	return F
end

# ╔═╡ 7689e15c-a631-431e-89b2-a5833e434fc4
"Generates the matrix to be considered in the power method"
function gen_Lq(q::Real)
	return stack(2*real(gen_approx_k0(q,1) + gen_approx_F(q,2)),dims=1)
end

# ╔═╡ 69eb6d7a-d477-4fa3-825b-70e032388149
"Uses the power method to return an approximation of the largest eigenvalue"
function power_method(Lq)
	v = ones((Nc+1)^2)
	for _ in 0:500
		v = (1/norm(v))*(Lq*v)
	end
	return dot(v,Lq*v)/dot(v,v)
end

# ╔═╡ bb02445d-34a1-4344-a1fe-5f09cfbb7f0d
"Implements the secant method. kinda https://mmas.github.io/secant-julia"
function secant()
	q1 = 1.3
	q2 = 1.31
	maxiter = 30
	tol = 0.001
    for _ in 1:maxiter
		h1 = power_method(gen_Lq(q1)) - 1
		h2 = power_method(gen_Lq(q2)) - 1
        q = q1 - (h1*(q2-q1))/(h2-h1)
        if abs(q-q2) < tol
            return q
        end
        q1 = q2
        q2 = q
    end
    error("Max iteration exceeded")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[compat]
SpecialFunctions = "~2.3.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.4"
manifest_format = "2.0"
project_hash = "dfc7de891366ddc30d160a79ab25207f28ae8702"

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

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

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

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

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

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

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
# ╠═236bdd31-4496-48e0-be58-7319738fc1b6
# ╠═6be3c15b-e874-48b2-9d98-e7b8c2f2621a
# ╠═9144b032-ed27-4910-bf0e-b26f1bf3b8ff
# ╠═6a4be131-7112-4a03-b30a-e6d1fe43e0f1
# ╠═1041df01-3192-4f28-a0d9-77733cbb6139
# ╠═aa71d4ec-94aa-4b43-a6cb-1361948d9dc7
# ╠═6a94dbda-1ee4-484c-a7cf-2c58aba99ef5
# ╠═84fcb6c5-dcb7-49d4-95d8-33aabe5a5a55
# ╠═5afdaf31-2bfe-4b1f-8991-a36e1dba5811
# ╠═c8f312cf-b086-4307-8a06-7ff3ef762a56
# ╠═77576ece-4e29-46e5-8dc2-35d87ee74449
# ╠═0c77f33d-3d6b-4fba-b80e-a5a85460c973
# ╠═626dbef0-bc4c-49f1-9eb0-fe9a08fac272
# ╠═746ea648-f98d-47b9-99e8-9655d55186f7
# ╠═c03fd4c5-ad48-4194-9015-ef44b556985f
# ╠═f1bcffef-817d-44e2-ab36-fc95dfe867e1
# ╠═e66b214d-d997-4afd-b140-062b7938ad7e
# ╠═609b7838-dd9d-4668-af75-25c85d3c0ced
# ╠═7689e15c-a631-431e-89b2-a5833e434fc4
# ╠═69eb6d7a-d477-4fa3-825b-70e032388149
# ╠═bb02445d-34a1-4344-a1fe-5f09cfbb7f0d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
