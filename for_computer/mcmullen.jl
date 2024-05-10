### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils
using Printf

# ╔═╡ e3e47920-d41f-11ee-0589-27b9f5e84d19
using LinearAlgebra

# ╔═╡ e0db369a-145b-49a5-88f4-92e2a4dfc0bc
using SparseArrays

"transforms a circle configuration with two parallel horizontal lines into a configuration inside the unit circle"
function rescale(circles,circles2 = circles)
	s=1
	for v in circles2
		if v[2]==0&&v[3]==0&&v[4]==1
			s=v[1]/2
		end
	end
	return [[[0.5,.5,-0.5,0] [.5,0.5,.5,0] [0,0,0,1] [1,-1,0,0]]*[[1/s,0,0,0] [0,s,0,0] [0,0,1,0] [0,0,0,1]]*v for v in circles]
end

# ╔═╡ 0a0798b9-2e24-481f-91d0-a7d513936c73
"transforms a circle configuration with two parallel vertical lines into a configuration with two parallel horizontal lines"
function reflect_xy(circles)
	return [[[1,0,0,0] [0,1,0,0] [0,0,0,1] [0,0,1,0]]*v for v in circles]
end

# ╔═╡ e6546fce-adb6-47f0-8a1c-a144e49cd10e
"This is the data tabulated by Chait, Cui, and Stier, for polyhedra with at most seven vertices and at most seven edges."

# polyhedron coords = rescale(polydual,poly)
# dual polyhedron coords = rescale(reflectxy(poly),reflectxy(polydual))

# ╔═╡ 6c12c61b-f87f-4143-8b12-e94935469aea
tetrahedron=rescale([[4,1,1,2],[0,1,1,0],[4,0,1,0],[0,0,-1,0]],[[4,1,2,1],[0,1,0,1],[4,0,0,1],[0,0,0,-1]])

# ╔═╡ b750135d-0d45-4129-abc0-a43f0fde2221
square_pyramid=rescale([[4,2,2sqrt(2),1],[8,1,2sqrt(2),1],[4,0,0,1],[0,1,0,1],[0,0,0,-1]])

# ╔═╡ 90b4d222-6a1e-4961-8011-c570c33b8acc
triangular_bipyramid=rescale([[64,1,4,7],[16,1/4,2,1],[0,1/4,0,1],[16,0,0,1],[0,0,0,-1]])

# ╔═╡ 1c1d0921-d0d3-4b89-9878-bea634c6d6dd
triangular_prism=rescale(reflect_xy([[32,1/2,1,4],[48,1/2,3,4],[0,1/4,1,0],[32,3/4,3,4],[0,0,-1,0],[16,0,1,0]]))

# ╔═╡ a43e7624-f98e-4677-a83d-fe8ecd394851
sixv6f1=rescale([[2(-1+sqrt(5)),1,sqrt(2(-1+sqrt(5))),1],[6-2sqrt(5),2+sqrt(5),sqrt(2(1+sqrt(5))),1],[0,1/2*(3+sqrt(5)),0,1],[2(-1+sqrt(5)),1/2*(1+sqrt(5)),0,sqrt(5)],[4,0,0,1],[0,0,0,-1]])

# ╔═╡ 000627cd-1022-4f51-a64b-94c87a02e6d8
pentagonal_pyramid=rescale([[0,1,0,1],[4,1/2*(3+sqrt(5)),1+sqrt(5),1],[6+2sqrt(5),1/2*(3+sqrt(5)),3+sqrt(5),1],[6+2sqrt(5),1,1+sqrt(5),1],[4,0,0,1],[0,0,0,-1]])

# ╔═╡ 9bc9ee4d-727d-4f15-9a35-17b77102437b
octahedron = rescale([[16*sqrt(2),sqrt(2),5,2*sqrt(2)],[8*sqrt(2),1/sqrt(2),1,2*sqrt(2)],[8*sqrt(2),sqrt(2),3,2*sqrt(2)],[0,1/sqrt(2),1,0],[8*sqrt(2),1/sqrt(2),3,0],[16*sqrt(2),1/sqrt(2),3,2*sqrt(2)],[0,0,-1,0],[8*sqrt(2),0,1,0]],[[16,1,2*sqrt(2),3],[16,1/2,2*sqrt(2),1],[8,1,2*sqrt(2),1],[0,1/2,0,1],[8,0,0,1],[0,0,0,-1]])

# ╔═╡ 0d93c5a9-2805-4bf8-b59c-78c4beb455a7
cube = rescale(reflect_xy([[16,1,2*sqrt(2),3],[16,1/2,2*sqrt(2),1],[8,1,2*sqrt(2),1],[0,1/2,0,1],[8,0,0,1],[0,0,0,-1]]),reflect_xy([[16*sqrt(2),sqrt(2),5,2*sqrt(2)],[8*sqrt(2),1/sqrt(2),1,2*sqrt(2)],[8*sqrt(2),sqrt(2),3,2*sqrt(2)],[0,1/sqrt(2),1,0],[8*sqrt(2),1/sqrt(2),3,0],[16*sqrt(2),1/sqrt(2),3,2*sqrt(2)],[0,0,-1,0],[8*sqrt(2),0,1,0]]))

# ╔═╡ e7a26e3e-bf9f-4c82-aded-d6086f8a5996
sixv7f1=rescale([[4sqrt(2),1+sqrt(2),1+2sqrt(2),0],[4sqrt(2),3+2sqrt(2),3+2sqrt(2),0],[0,1+sqrt(2),1,0],[4,1,1,2],[2sqrt(2),1+sqrt(2),1+sqrt(2),sqrt(2)],[0,0,-1,0],[4,0,1,0]],[[4,1,2,1],[8,3+2sqrt(2),2(2+sqrt(2)),1],[4,3+2sqrt(2),2(1+sqrt(2)),1],[0,1,0,1],[4,0,0,1],[0,0,0,-1]])

# ╔═╡ d78ce5a0-dfb2-47d9-83f4-405a4e255124
sixv7f1_dual=rescale(reflect_xy([[4,1,2,1],[8,3+2sqrt(2),2(2+sqrt(2)),1],[4,3+2sqrt(2),2(1+sqrt(2)),1],[0,1,0,1],[4,0,0,1],[0,0,0,-1]]),reflect_xy([[4sqrt(2),1+sqrt(2),1+2sqrt(2),0],[4sqrt(2),3+2sqrt(2),3+2sqrt(2),0],[0,1+sqrt(2),1,0],[4,1,1,2],[2sqrt(2),1+sqrt(2),1+sqrt(2),sqrt(2)],[0,0,-1,0],[4,0,1,0]]))

# ╔═╡ 3e736c46-94fc-4578-8103-aae448e22431
sixv7f2=rescale([[0,sqrt(3),1,0],[4sqrt(3),2sqrt(3),5,0],[20/sqrt(3),8/sqrt(3),7,4/sqrt(3)],[8/sqrt(3),2/sqrt(3),1,4/sqrt(3)],[8/sqrt(3),5/sqrt(3),3,4/sqrt(3)],[0,0,-1,0],[2sqrt(3),0,1,0]],[[0,1,0,1],[4,3,2sqrt(3),1],[12,4,4sqrt(3),1],[16/3,4/3,4/sqrt(3),5/3],[4,0,0,1],[0,0,0,-1]])

# ╔═╡ 6753010e-e6cf-48aa-99d3-cbc65df1c1c7
sixv7f2_dual=rescale(reflect_xy([[0,1,0,1],[4,3,2sqrt(3),1],[12,4,4sqrt(3),1],[16/3,4/3,4/sqrt(3),5/3],[4,0,0,1],[0,0,0,-1]]),reflect_xy([[0,sqrt(3),1,0],[4sqrt(3),2sqrt(3),5,0],[20/sqrt(3),8/sqrt(3),7,4/sqrt(3)],[8/sqrt(3),2/sqrt(3),1,4/sqrt(3)],[8/sqrt(3),5/sqrt(3),3,4/sqrt(3)],[0,0,-1,0],[2sqrt(3),0,1,0]]))

# ╔═╡ f40f626c-a421-4364-a20f-b76c28c4a287
octahedron_orig = [[16,1,2sqrt(2),3],[16,1/2,2sqrt(2),1],[8,1,2sqrt(2),1],[0,1/2,0,1],[8,0,0,1],[0,0,0,-1]]

# ╔═╡ 230e3a15-6d0f-44fc-b48f-4d16dbce70d3
sevenv7f1=rescale([[0,1.306562965,1.000000000,0],[5.226251860,2.230442497,2.414213562,2.613125930],[7.391036260,1.847759065,1.000000000,3.695518130],[20.00832438,2.230442497,2.414213562,6.308644060],[20.90500744,1.306562965,1.000000000,5.226251860],[0,0,-1.000000000,0],[3.061467459,0,1.000000000,0]],[[8,1,0,3],[0,1,0,1],[4,1+1/sqrt(2),sqrt(2(2+sqrt(2))),1],[2(2+sqrt(2)),1/2+1/sqrt(2),sqrt(2+sqrt(2)),1+sqrt(2)],[12+8sqrt(2),1+1/sqrt(2),sqrt(2(2+sqrt(2))),3+2sqrt(2)],[8,0,0,1],[0,0,0,-1]])

# ╔═╡ 7f4ad4d2-6138-4d8d-b883-0ba8f9755752
sevenv7f1_dual=rescale(reflect_xy([[8,1,0,3],[0,1,0,1],[4,1+1/sqrt(2),sqrt(2(2+sqrt(2))),1],[2(2+sqrt(2)),1/2+1/sqrt(2),sqrt(2+sqrt(2)),1+sqrt(2)],[12+8sqrt(2),1+1/sqrt(2),sqrt(2(2+sqrt(2))),3+2sqrt(2)],[8,0,0,1],[0,0,0,-1]]),reflect_xy([[0,1.306562965,1.000000000,0],[5.226251860,2.230442497,2.414213562,2.613125930],[7.391036260,1.847759065,1.000000000,3.695518130],[20.00832438,2.230442497,2.414213562,6.308644060],[20.90500744,1.306562965,1.000000000,5.226251860],[0,0,-1.000000000,0],[3.061467459,0,1.000000000,0]]))

# ╔═╡ 3d3592ba-58c7-43f4-ab62-84cd8b966120
sevenv7f2=rescale([[8.65914,0.270598,1.,1.53073],[20.905,0.653281,3.82843,0.],[71.3742,3.53701,15.4853,3.69552],[56.5921,2.23044,10.6569,3.69552],[0.,0.46194,1.,0.],[0.,0.,-1.,0.],[14.7821,0.,1.,0.]],[[16,1,2sqrt(2+sqrt(2)),-1+2sqrt(2)],[0,1/(2sqrt(2)),0,1],[16(1+sqrt(2)),1+3/(2sqrt(2)),2sqrt(10+7sqrt(2)),1],[24(2+sqrt(2)),2+sqrt(2),1/382*(2035+sqrt(17610545)),1+2sqrt(2)],[8(1+sqrt(2)),1/(2sqrt(2)),sqrt(2(2+sqrt(2))),1],[8sqrt(2),0,0,1],[0,0,0,-1]])

# ╔═╡ 4cbf5fa9-9d99-4831-a622-e230eedbfb9d
sevenv7f2_dual=rescale(reflect_xy([[16,1,2sqrt(2+sqrt(2)),-1+2sqrt(2)],[0,1/(2sqrt(2)),0,1],[16(1+sqrt(2)),1+3/(2sqrt(2)),2sqrt(10+7sqrt(2)),1],[24(2+sqrt(2)),2+sqrt(2),1/382*(2035+sqrt(17610545)),1+2sqrt(2)],[8(1+sqrt(2)),1/(2sqrt(2)),sqrt(2(2+sqrt(2))),1],[8sqrt(2),0,0,1],[0,0,0,-1]]),reflect_xy([[8.65914,0.270598,1.,1.53073],[20.905,0.653281,3.82843,0.],[71.3742,3.53701,15.4853,3.69552],[56.5921,2.23044,10.6569,3.69552],[0.,0.46194,1.,0.],[0.,0.,-1.,0.],[14.7821,0.,1.,0.]]))

# ╔═╡ 8c20ba94-2576-44cd-916b-e9f2c83148c0
sevenv7f3=rescale([[14+6sqrt(5),1,2+sqrt(5),1+sqrt(5)],[4(2+sqrt(5)),-2+sqrt(5),1,2],[0,(3-sqrt(5))/2,1,0],[14+6sqrt(5),(-1+sqrt(5))/2,2+sqrt(5),0],[8(2+sqrt(5)),(-1+sqrt(5))/2,2+sqrt(5),2],[0,0,-1,0],[4(2+sqrt(5)),0,1,0]],[[0,1,0,1],[4,3,2sqrt(3),1],[12,4,4sqrt(3),1],[16/3,4/3,4/sqrt(3),5/3],[4,0,0,1],[0,0,0,-1]])

# ╔═╡ c2544c12-fe91-44d8-a280-0cdcfc52a98b
sevenv7f3_dual=rescale(reflect_xy([[0,1,0,1],[4,3,2sqrt(3),1],[12,4,4sqrt(3),1],[16/3,4/3,4/sqrt(3),5/3],[4,0,0,1],[0,0,0,-1]]),reflect_xy([[14+6sqrt(5),1,2+sqrt(5),1+sqrt(5)],[4(2+sqrt(5)),-2+sqrt(5),1,2],[0,(3-sqrt(5))/2,1,0],[14+6sqrt(5),(-1+sqrt(5))/2,2+sqrt(5),0],[8(2+sqrt(5)),(-1+sqrt(5))/2,2+sqrt(5),2],[0,0,-1,0],[4(2+sqrt(5)),0,1,0]]))

# ╔═╡ 84b36d9f-6dbb-42bb-9262-b4b012112a12
sevenv7f4=rescale([[128,1/2,1,8],[144,1/2,3,8],[0,1/4,1,0],[16,1/4,1,2],[144,3/4,3,10],[16,0,1,0],[0,0,-1,0]],[[240,1,4,15],[48,1/4,2,3],[16,1/4,2,1],[0,1/4,0,1],[32,1/4,0,3],[32,0,0,1],[0,0,0,-1]])

# ╔═╡ 6c51b8f1-130a-4920-948b-c022f59a20b5
sevenv7f4_dual=rescale(reflect_xy([[240,1,4,15],[48,1/4,2,3],[16,1/4,2,1],[0,1/4,0,1],[32,1/4,0,3],[32,0,0,1],[0,0,0,-1]]),reflect_xy([[128,1/2,1,8],[144,1/2,3,8],[0,1/4,1,0],[16,1/4,1,2],[144,3/4,3,10],[16,0,1,0],[0,0,-1,0]]))

# ╔═╡ c0a7f290-8a3e-4317-ae21-01f0ffdd4d54
sevenv7f5=rescale([[261.837,0.968135,10.0696,12.3732],[54.4322,0.174917,1.,3.08563],[0.,0.174917,1.,0.],[143.138,0.726701,6.25932,8.11415],[54.4322,0.416351,3.76056,3.08563],[0.,0.,-1.,0.],[22.868,0.,1.,0.]],[[243.972,1.,8.96427,12.8302],[0.,0.113375,0.,1.],[35.2811,0.269865,3.08563,1.],[83.9787,0.372488,4.25903,3.76056],[146.577,0.411513,4.70524,6.25932],[35.2811,0.,0.,1.],[0.,0.,0.,-1.]])

# ╔═╡ b143e1ed-f708-44a1-9cf4-12e5ab764ccd
sevenv7f5_dual=rescale(reflect_xy([[243.972,1.,8.96427,12.8302],[0.,0.113375,0.,1.],[35.2811,0.269865,3.08563,1.],[83.9787,0.372488,4.25903,3.76056],[146.577,0.411513,4.70524,6.25932],[35.2811,0.,0.,1.],[0.,0.,0.,-1.]]),reflect_xy([[261.837,0.968135,10.0696,12.3732],[54.4322,0.174917,1.,3.08563],[0.,0.174917,1.,0.],[143.138,0.726701,6.25932,8.11415],[54.4322,0.416351,3.76056,3.08563],[0.,0.,-1.,0.],[22.868,0.,1.,0.]]))

# ╔═╡ 4f99384f-c825-4334-b161-925bf0faaa40
sevenv7f6=rescale([[4.60386,3.82663,3.64944,2.30193],[3.47535,0.868837,1.,1.73767],[16.8015,3.5445,7.15919,3.04941],[3.47535,1.5247,2.50976,0.],[0.,2.0198,1.,0.],[0.,0.,-1.,0.],[4.60386,0.,1.,0]],[[0.,1.,0.,1.],[4.,4.0796,4.0396,1.],[5.29887,2.0796,3.04941,1.64944],[12.3184,2.32472,5.35133,1.],[12.3184,1.75488,4.0396,2.50976],[4.,0.,0.,1.],[0.,0.,0.,-1.]])

# ╔═╡ 1d614f85-32d3-473d-b436-2f5f0fc8ce02
sevenv7f6_dual=rescale(reflect_xy([[0.,1.,0.,1.],[4.,4.0796,4.0396,1.],[5.29887,2.0796,3.04941,1.64944],[12.3184,2.32472,5.35133,1.],[12.3184,1.75488,4.0396,2.50976],[4.,0.,0.,1.],[0.,0.,0.,-1.]]),reflect_xy([[4.60386,3.82663,3.64944,2.30193],[3.47535,0.868837,1.,1.73767],[16.8015,3.5445,7.15919,3.04941],[3.47535,1.5247,2.50976,0.],[0.,2.0198,1.,0.],[0.,0.,-1.,0.],[4.60386,0.,1.,0]]))

# ╔═╡ bef41f35-3c1e-430a-8e41-61075a736bc9
sevenv7f7=rescale([[53.5306,0.255122,1.,3.69552],[69.2094,0.690643,3.,6.30864],[15.6788,0.435521,1.,2.61313],[0.,0.255122,1.,0.],[129.234,0.871042,5.82843,8.92177],[0.,0.,-1.,0.],[15.6788,0.,1.,0.]],[[72+48sqrt(2),1,2sqrt(2(2+sqrt(2))),5+4sqrt(2)],[6(2+sqrt(2)),1/(3sqrt(2)),0,1+sqrt(2)],[0,1/3,0,1],[12,(2-sqrt(2))/3,2sqrt(2-sqrt(2)),1],[48+36sqrt(2),sqrt(2)/3,2sqrt(2+sqrt(2)),3+2sqrt(2)],[12(1+sqrt(2)),0,0,1],[0,0,0,-1]])

# ╔═╡ f612362e-83a3-44fd-a9cf-a2816308c500
sevenv7f7_dual=rescale(reflect_xy([[72+48sqrt(2),1,2sqrt(2(2+sqrt(2))),5+4sqrt(2)],[6(2+sqrt(2)),1/(3sqrt(2)),0,1+sqrt(2)],[0,1/3,0,1],[12,(2-sqrt(2))/3,2sqrt(2-sqrt(2)),1],[48+36sqrt(2),sqrt(2)/3,2sqrt(2+sqrt(2)),3+2sqrt(2)],[12(1+sqrt(2)),0,0,1],[0,0,0,-1]]),reflect_xy([[53.5306,0.255122,1.,3.69552],[69.2094,0.690643,3.,6.30864],[15.6788,0.435521,1.,2.61313],[0.,0.255122,1.,0.],[129.234,0.871042,5.82843,8.92177],[0.,0.,-1.,0.],[15.6788,0.,1.,0.]]))


# ╔═╡ 2e9b6f4e-025f-4945-9a00-4e2ab8aea517
"For reformatting data copied from Mathematica"
function convert_to_matrix(string)
	str="[["
	for i in 1:length(string)
		if string[i]==' '
		elseif string[i]=='['
			str=str*"("
		elseif string[i]==']'
			str=str*")"
		elseif string[i]=='S'
			str=str*"s"
		elseif string[i]=='\t'
			str=str*","
		elseif string[i]=='\n'
			str=str*"],["
		elseif string[i]=='{'
			str=str*"["
		elseif string[i]=='}'
			str=str*"]"
		else
			str=str*string[i]
		end
	end
	str=str*"]]"
	return str
end

"reflect circle2 over circle1"
function reflect(c1, c2)
	return c2 - c1*(transpose(c1)*[[0, -1, 0, 0] [-1, 0, 0, 0] [0, 0, 2, 0] [0, 0, 0, 2]]*c2)
end

tetra = [[0, 0, 0, -1],[1, 1, 1, 1],[1, 1, -1, 1],[4, 0, 0, 1]]

struct Dual
	coords
	word
	children
	parent
end

numbers = Dict()

whole_tree = Dict()

function find_siblings(node)
	par = node.parent
	sibs = [par]
	for child in par.children
		if child.coords != node.coords
			push!(sibs,child)
		end
	end
	return sibs
end

function find_leaves(node)
	current = [node]
	new = []
	while length(current) > 0
		leaf = pop!(current)
		if length(leaf.children) == 0
			push!(new,leaf)
		else 
			for child in leaf.children
				push!(current, child)
			end
		end
	end
	return new
end

function spawn(circ)
	sibs = find_siblings(circ)
	child_word = push!(copy(circ.word),last(sibs[1].word))
	child_coords = -1*reflect(circ.coords,sibs[1].coords)
	child = Dual(child_coords,child_word,[],circ)
	push!(circ.children,child)
	whole_tree[join(child_word)] = child
	for sib in sibs[2:end]
		child_word = push!(copy(circ.word),last(sib.word))
		child_coords = reflect(circ.coords,sib.coords)
		child = Dual(child_coords,child_word,[],circ)
		push!(circ.children,child)
		whole_tree[join(child_word)] = child
	end
	return circ
end

function sample_value(leaf,func)
	point = leaf.coords[3]/leaf.coords[1] + 1im*leaf.coords[4]/leaf.coords[1]
	if dual_circles[func][1] == 0
		return 1
	else
		center = dual_circles[func][3]/dual_circles[func][1] + 1im*dual_circles[func][4]/dual_circles[func][1]
		return abs(((point-center)*dual_circles[func][1])^(-2))
	end
end

function generate_tree()
	stem = []
	tree = [] 
	for i in 1:N
		start = Dual(dual_circles[i],[i],[],0)
		push!(stem,start)
		whole_tree[join([i])] = start
	end
	for (i,root) in enumerate(stem)
		for j in 1:N
			if j != i
				child_word = [i,j]
				child = Dual(reflect(root.coords,dual_circles[j]),child_word,[],root)
				push!(root.children,child)
				push!(tree,child)
				whole_tree[join(child_word)] = child
			end
		end
	end
	leaves = []
	while length(tree) > 0
		node = pop!(tree)
		if node.coords[1] > ceiling 
			push!(leaves,node)
		else
			node = spawn(node)
			for child in node.children
				push!(tree,child) 
			end
		end
	end
	
	for (i, leaf) in enumerate(leaves)
		numbers[join(leaf.word)] = i
	end

	return leaves
end

function generate_data()
	col = []
	row = []
	nonzero = []
	for (i,node) in enumerate(data)
		sample = whole_tree[join(node.word[2:end])]
		for leaf in find_leaves(sample)
			push!(row, i)
			push!(col, numbers[join(leaf.word)])
			push!(nonzero, sample_value(leaf,first(node.word)))
		end
	end
	return col, row, nonzero, length(data)
end

function generate_data_simplified()
	col = []
	row = []
	nonzero = []
	for (i,node) in enumerate(data)
		sample = whole_tree[join(node.word[2:end])]
		for leaf in find_leaves(sample)
			push!(row, i)
			push!(col, numbers[join(leaf.word)])
			push!(nonzero, BigFloat(leaf.coords[1]/reflect(dual_circles[node.word[1]],leaf.coords)[1]))
		end
	end
	return col, row, nonzero, length(data)
end

function power_method(col,row,data,l,a)
	matr = sparse(col, row, data.^a)
    vec = ones(BigFloat, l)
    previous_entry = vec[1]
    previous_val = 0
    current = matr * vec
    current_val = current[1] / previous_entry
    count = 0
    while count < 1000 && abs(current_val - previous_val) > 1e-8
        previous_val = current_val
        previous_entry = current[1]
        current = matr * current
        current_val = current[1] / previous_entry
        count += 1
	end
    @printf("power method iterations: %d\n", count)
    return current_val
end

function secant(x0,y0,x1,y1,z)
	return x0 - (y0-z) * ((x1-x0)/(y1-y0))
end


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

ceiling = parse(Float64, ARGS[1])

dual_circles = tetra

N = size(dual_circles)[1]

t1 = time()

data = generate_tree()

col2,row2,nonzero2,l2 = generate_data()
@printf("ceiling: %d\n", ceiling)

@printf("partitions: %d\n", l2)

n,q,v = secant_method(col2,row2,nonzero2,l2,1,1.3,1.31,1e-8,100)

@printf("dimension: %f, eigenval: %f, iterations: %d\n", q, v, n)

@printf("time: %f\n", (time()-t1)/60)

@printf("\n")
