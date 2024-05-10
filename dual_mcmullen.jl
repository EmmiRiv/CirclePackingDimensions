### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ e3e47920-d41f-11ee-0589-27b9f5e84d19
using LinearAlgebra

# ╔═╡ 21cf04e2-fdf9-4d4f-b1dd-c92a83594105
using Luxor

# ╔═╡ e0db369a-145b-49a5-88f4-92e2a4dfc0bc
using SparseArrays

# ╔═╡ 1160700e-81fe-4161-aaa0-4362fcf69e49
using Printf

# ╔═╡ d6be8af4-aed1-4897-811c-722e610dc0f6
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

# ╔═╡ 9434ac79-95bc-4f41-8caa-9a7cce05ce67
hexagonal_pyramid=rescale([[12,1,2sqrt(3),1],[16,3,4sqrt(3),1],[12,4,4sqrt(3),1],[4,3,2sqrt(3),1],[0,1,0,1],[4,0,0,1],[0,0,0,-1]])

# ╔═╡ f0fb32f7-5434-4f11-bc60-3e311a03fc86
"Here is an illustration of each circle configuration."

# ╔═╡ 8965d128-b902-4d10-aa57-da0f896f8e81
function label(v)
	if v[1]>0
		text(string(v[1]), 300*v[3]/abs(v[1]), 300*v[4]/abs(v[1]), halign=:center, valign=:center)
	end
end

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

# ╔═╡ 59639686-c6b4-4304-8c24-e203c14d3fb0
tetra = [[0, 0, 0, -1],[1, 1, 1, 1],[1, 1, -1, 1],[4, 0, 0, 1]]

# ╔═╡ 14bb4ee3-7059-44b6-98dc-f045ba1dd76a
dual_circles = cube

# ╔═╡ 587ac08f-25a9-4d38-be4a-484a6cb0468c
N = size(dual_circles)[1]

# ╔═╡ eb8ac9dc-ebf5-474e-9534-f4301565f761
struct Dual
	coords
	word
	children
	parent
end

# ╔═╡ 37adfcd7-b9c6-417f-a146-15f9db7abda1
numbers = Dict()

# ╔═╡ 477f7ecf-d428-4dd3-bf56-45c4bb834ecc
whole_tree = Dict()

# ╔═╡ 18ceb249-cdb9-466d-ab15-7ce50ba3c143
"reflect circle2 over circle1"
function reflect(c1, c2)
	return c2 - c1*(transpose(c1)*[[0, -1, 0, 0] [-1, 0, 0, 0] [0, 0, 2, 0] [0, 0, 0, 2]]*c2)
end

# ╔═╡ c4e483ae-defd-4951-bc02-27994518b273
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

# ╔═╡ 6e090376-4512-48fe-8387-fed03fe7ca10
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

# ╔═╡ c7f09a39-a844-4dcf-8bb0-404438fa8be6
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

# ╔═╡ a4411349-54d1-4e91-8d0d-2c0a7c18932c
function sample_value(leaf,func)
	point = leaf.coords[3]/leaf.coords[1] + 1im*leaf.coords[4]/leaf.coords[1]
	if dual_circles[func][1] == 0
		return 1
	else
		center = dual_circles[func][3]/dual_circles[func][1] + 1im*dual_circles[func][4]/dual_circles[func][1]
		return abs(((point-center)*dual_circles[func][1])^(-2))
	end
end

# ╔═╡ 456f05a9-f100-4804-8350-483dc9e249e4
function power_method(col,row,data,l,a)
	matr = sparse(col, row, data.^a)
    vec = ones(BigFloat, l)
    previous_entry = vec[1]
    previous_val = 0
    current = matr * vec
    current_val = current[1] / previous_entry
    count = 0
    while count < 1000 && abs(current_val - previous_val) > 1e-6
        previous_val = current_val
        previous_entry = current[1]
        current = matr * current
        current_val = current[1] / previous_entry
        count += 1
	end
    return current_val
end

# ╔═╡ 4156ae68-7cc0-44e5-a061-4e386a84b82d
ceiling = 20

# ╔═╡ 17554d13-9fb3-42ca-8c8a-7845121ab7d0
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
	while length(tree)>0
		
		node = pop!(tree)
		if node.coords[1]<0
			return node
		end
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

# ╔═╡ 510c0487-3843-4542-9732-8018c6bd4f85
data = generate_tree()

# ╔═╡ 3474f2b9-2643-435a-8a14-8386eca94864
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

# ╔═╡ 87a02f27-99f3-44c9-b59c-cc5743522ab7
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

# ╔═╡ ccda2545-bcd7-4a08-80f0-4deaccbdc007
find_siblings(data.parent)

# ╔═╡ c9f2706f-8781-49da-8ee7-56337818e7cc
data.parent

# ╔═╡ 40b9e184-a409-4872-905a-7b115767fa37
# ╠═╡ disabled = true
#=╠═╡
col,row,nonzero,l = generate_data_simplified()
  ╠═╡ =#

# ╔═╡ e185c3ee-8e08-47a9-93ab-120dd7101c42
# ╠═╡ disabled = true
#=╠═╡
col2,row2,nonzero2,l2 = generate_data()
  ╠═╡ =#

# ╔═╡ 4524df63-5985-46ec-85d8-c364668a7cb1
# ╠═╡ disabled = true
#=╠═╡
power_method(col,row,nonzero,l,BigFloat(1.305686728))
  ╠═╡ =#

# ╔═╡ b4cb2d45-ba17-4846-aefc-0bf7033529f6
# ╠═╡ disabled = true
#=╠═╡
power_method(col2,row2,nonzero2,l2,BigFloat(1.305686728))
  ╠═╡ =#

# ╔═╡ 5b99bdc5-4010-4f20-9ecd-ba7af405176c
# ╠═╡ disabled = true
#=╠═╡
secant_method(col,row,nonzero,l,1,3.0,3.1,1e-4,10)
  ╠═╡ =#

# ╔═╡ c97e6297-e53f-4437-984f-f2c67ec2cc12
# ╠═╡ disabled = true
#=╠═╡
secant_method(col2,row2,nonzero2,l2,1,3.0,3.1,1e-6,100)
  ╠═╡ =#

# ╔═╡ c47e0b49-0748-45d3-90f5-bd02d05a45ac
function secant(x0,y0,x1,y1,z)
	return x0 - (y0-z) * ((x1-x0)/(y1-y0))
end

# ╔═╡ a0b1115e-d614-4595-84fe-8fb4de91ec39
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

# ╔═╡ 25642183-fa77-43b4-adb5-7fc5346f6ddc
# ╠═╡ disabled = true
#=╠═╡
@drawsvg begin
	# for ceiling = 50
	background("white")
	sethue("black")
	s = 300
	for c in dual_circles
		draw_circle(c)
	end

	thing = "4143"
	
	sethue("olivedrab")
	c = whole_tree[thing].coords
	draw_circle(c)

	sethue("chartreuse4")
	c = whole_tree[thing[2:end]].coords
	draw_circle(c)
	
	for circ in find_leaves(whole_tree[thing[2:end]])
		sethue("chartreuse3")
		d = circ.coords
		circle(Point(s*d[3]/d[1],s*d[4]/d[1]),s/d[1],:stroke)
		sethue("steelblue1")
		circle(Point(s*d[3]/d[1],s*d[4]/d[1]),1,:fill)
		v = reflect(dual_circles[4],d)
		sethue("steelblue")
		circle(Point(s*v[3]/v[1],s*v[4]/v[1]),0.75,:fill)
	end
		
end
  ╠═╡ =#

# ╔═╡ 76f27c49-71e0-494e-b3a8-bfac6821c23b
function draw_circle(v)
	sc = 50
	if v[1]!=0
		circle(2*sc*v[3]/v[1], 2*sc*v[4]/v[1], 2*sc/abs(v[1]), :stroke)
	else
		rule(Point(sc*v[2]*v[3], sc*v[2]*v[4]), atan(v[3],-v[4]))
	end
end

# ╔═╡ c59f07cb-07ea-4e48-ba7c-f40297f1e5e8
@drawsvg begin
	background("white")
	sethue("blue")
	for c in tetrahedron
		draw_circle(c)
	end
end

# ╔═╡ 0044c7b7-60b3-4031-9a08-5fa53aba27fd
@svg begin
	background("white")
	for thing in data
		setcolor("black")
		draw_circle(thing.coords)
	end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Luxor = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
Luxor = "~3.8.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.4"
manifest_format = "2.0"
project_hash = "93a3d29a3e2f0c3acf5c8dd97301dea30f5feb16"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "886826d76ea9e72b35fcd000e535588f7b60f21d"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60b1194df0a3298f460063de985eae7b01bc011a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.1+0"

[[deps.Juno]]
deps = ["Base64", "Logging", "Media", "Profile"]
git-tree-sha1 = "07cb43290a840908a771552911a6274bc6c072c7"
uuid = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
version = "0.8.4"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pango_jll", "Pkg", "gdk_pixbuf_jll"]
git-tree-sha1 = "ae0923dab7324e6bc980834f709c4cd83dd797ed"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.54.5+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Luxor]]
deps = ["Base64", "Cairo", "Colors", "DataStructures", "Dates", "FFMPEG", "FileIO", "Juno", "LaTeXStrings", "PrecompileTools", "Random", "Requires", "Rsvg"]
git-tree-sha1 = "aa3eb624552373a6204c19b00e95ce62ea932d32"
uuid = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc"
version = "3.8.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Media]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "75a54abd10709c01f1b86b84ec225d26e840ed58"
uuid = "e89f7d12-3494-54d1-8411-f7d8b9ae1f27"
version = "0.5.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc6e1927ac521b659af340e0ca45828a3ffc748f"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.12+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4745216e94f71cb768d58330b059c9b76f32cb66"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.14+0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rsvg]]
deps = ["Cairo", "Glib_jll", "Librsvg_jll"]
git-tree-sha1 = "3d3dc66eb46568fb3a5259034bfc752a0eb0c686"
uuid = "c4c386cf-5103-5370-be45-f3a111cca3b8"
version = "1.0.0"

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

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "801cbe47eae69adc50f36c3caec4758d2650741b"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.2+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "e9190f9fb03f9c3b15b9fb0c380b0d57a3c8ea39"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.8+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "93284c28274d9e75218a416c65ec49d0e0fcdf3d"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.40+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╠═e3e47920-d41f-11ee-0589-27b9f5e84d19
# ╠═21cf04e2-fdf9-4d4f-b1dd-c92a83594105
# ╠═e0db369a-145b-49a5-88f4-92e2a4dfc0bc
# ╠═1160700e-81fe-4161-aaa0-4362fcf69e49
# ╠═d6be8af4-aed1-4897-811c-722e610dc0f6
# ╠═0a0798b9-2e24-481f-91d0-a7d513936c73
# ╠═e6546fce-adb6-47f0-8a1c-a144e49cd10e
# ╠═6c12c61b-f87f-4143-8b12-e94935469aea
# ╠═c59f07cb-07ea-4e48-ba7c-f40297f1e5e8
# ╠═b750135d-0d45-4129-abc0-a43f0fde2221
# ╠═90b4d222-6a1e-4961-8011-c570c33b8acc
# ╠═1c1d0921-d0d3-4b89-9878-bea634c6d6dd
# ╠═a43e7624-f98e-4677-a83d-fe8ecd394851
# ╠═000627cd-1022-4f51-a64b-94c87a02e6d8
# ╠═9bc9ee4d-727d-4f15-9a35-17b77102437b
# ╠═0d93c5a9-2805-4bf8-b59c-78c4beb455a7
# ╠═e7a26e3e-bf9f-4c82-aded-d6086f8a5996
# ╠═d78ce5a0-dfb2-47d9-83f4-405a4e255124
# ╠═3e736c46-94fc-4578-8103-aae448e22431
# ╠═6753010e-e6cf-48aa-99d3-cbc65df1c1c7
# ╠═f40f626c-a421-4364-a20f-b76c28c4a287
# ╠═230e3a15-6d0f-44fc-b48f-4d16dbce70d3
# ╠═7f4ad4d2-6138-4d8d-b883-0ba8f9755752
# ╠═3d3592ba-58c7-43f4-ab62-84cd8b966120
# ╠═4cbf5fa9-9d99-4831-a622-e230eedbfb9d
# ╠═8c20ba94-2576-44cd-916b-e9f2c83148c0
# ╠═c2544c12-fe91-44d8-a280-0cdcfc52a98b
# ╠═84b36d9f-6dbb-42bb-9262-b4b012112a12
# ╠═6c51b8f1-130a-4920-948b-c022f59a20b5
# ╠═c0a7f290-8a3e-4317-ae21-01f0ffdd4d54
# ╠═b143e1ed-f708-44a1-9cf4-12e5ab764ccd
# ╠═4f99384f-c825-4334-b161-925bf0faaa40
# ╠═1d614f85-32d3-473d-b436-2f5f0fc8ce02
# ╠═bef41f35-3c1e-430a-8e41-61075a736bc9
# ╠═f612362e-83a3-44fd-a9cf-a2816308c500
# ╠═9434ac79-95bc-4f41-8caa-9a7cce05ce67
# ╠═f0fb32f7-5434-4f11-bc60-3e311a03fc86
# ╠═8965d128-b902-4d10-aa57-da0f896f8e81
# ╠═2e9b6f4e-025f-4945-9a00-4e2ab8aea517
# ╠═59639686-c6b4-4304-8c24-e203c14d3fb0
# ╠═14bb4ee3-7059-44b6-98dc-f045ba1dd76a
# ╠═587ac08f-25a9-4d38-be4a-484a6cb0468c
# ╠═eb8ac9dc-ebf5-474e-9534-f4301565f761
# ╠═37adfcd7-b9c6-417f-a146-15f9db7abda1
# ╠═477f7ecf-d428-4dd3-bf56-45c4bb834ecc
# ╠═18ceb249-cdb9-466d-ab15-7ce50ba3c143
# ╠═c4e483ae-defd-4951-bc02-27994518b273
# ╠═6e090376-4512-48fe-8387-fed03fe7ca10
# ╠═c7f09a39-a844-4dcf-8bb0-404438fa8be6
# ╠═17554d13-9fb3-42ca-8c8a-7845121ab7d0
# ╠═a4411349-54d1-4e91-8d0d-2c0a7c18932c
# ╠═3474f2b9-2643-435a-8a14-8386eca94864
# ╠═87a02f27-99f3-44c9-b59c-cc5743522ab7
# ╠═456f05a9-f100-4804-8350-483dc9e249e4
# ╠═4156ae68-7cc0-44e5-a061-4e386a84b82d
# ╠═0044c7b7-60b3-4031-9a08-5fa53aba27fd
# ╠═510c0487-3843-4542-9732-8018c6bd4f85
# ╠═ccda2545-bcd7-4a08-80f0-4deaccbdc007
# ╠═c9f2706f-8781-49da-8ee7-56337818e7cc
# ╠═40b9e184-a409-4872-905a-7b115767fa37
# ╠═e185c3ee-8e08-47a9-93ab-120dd7101c42
# ╠═4524df63-5985-46ec-85d8-c364668a7cb1
# ╠═b4cb2d45-ba17-4846-aefc-0bf7033529f6
# ╠═5b99bdc5-4010-4f20-9ecd-ba7af405176c
# ╠═c97e6297-e53f-4437-984f-f2c67ec2cc12
# ╠═c47e0b49-0748-45d3-90f5-bd02d05a45ac
# ╠═a0b1115e-d614-4595-84fe-8fb4de91ec39
# ╠═25642183-fa77-43b4-adb5-7fc5346f6ddc
# ╠═76f27c49-71e0-494e-b3a8-bfac6821c23b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
