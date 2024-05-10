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

dual_circles = octahedron

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
