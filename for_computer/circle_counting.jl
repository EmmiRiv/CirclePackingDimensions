### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 6a4e4863-05ad-401a-93df-ff757b6e3b6c
using LinearAlgebra

# ╔═╡ d3f5ce86-b727-42ba-b5b1-1194a29975a5
using Printf


"Fits a straight line through a set of points, `y = a₁ + a₂ * x`"
function linear_fit(x, y)

    
    sx = sum(x)
    sy = sum(y)

    m = length(x)

    sx2 = zero(sx.*sx)
    sy2 = zero(sy.*sy)
    sxy = zero(sx*sy)

    for i = 1:m
        sx2 += x[i]*x[i]
        sy2 += y[i]*y[i]
        sxy += x[i]*y[i]
    end

    a0 = (sx2*sy - sxy*sx) / ( m*sx2 - sx*sx )
    a1 = (m*sxy - sx*sy) / (m*sx2 - sx*sx)

    return (a0, a1)
end

# ╔═╡ 7a4cf9ec-c7b7-47af-b474-586aaaed97f0
inn(v, w)=-v[1]*w[2]-v[2]*w[1]+2*v[3]*w[3]+2*v[4]*w[4]

# ╔═╡ 722ca68c-8973-47b0-b54e-3119bf81489b
dual_circles=[[0, 0, 0, -1], [4, 0, 0, 1], [1,1,1,1], [1, 1,-1, 1]]

# ╔═╡ b4571059-8617-4d39-ab3e-909014fe2394
function all_circles(max)
	L=copy(dual_circles)
	i=1
	while i<=length(L)
		c=L[i]
		for d in dual_circles
			A=inn(c,d)
			if A<0&&c[1]-A*d[1]<max
				L=push!(L,c-A*d)
			end
		end
		i=i+1
	end
	return L
end

t1 = time()
ceiling = parse(Int64, ARGS[1])

@printf("ceiling: %d\n", ceiling)

circles = all_circles(ceiling)

hist=[0 for i in 1:ceiling+1]
for c in circles
    hist[abs(c[1])+1] += 1
end

tot = [0 for i in 0:ceiling+1]
for i in 1:ceiling+1
    tot[i+1]=tot[i]+hist[i]
end


logs = []
for i in 1:ceiling+1
    push!(logs,log(tot[i+1]))
end

dim = linear_fit([log(i) for i in div(ceiling,2):ceiling], [logs[i] for i in div(ceiling,2):ceiling])[2]
@printf("dimension: %f\n",dim)

@printf("circles counted: %d\n", trunc(ceiling^dim))

@printf("time: %f minutes\n", (time()-t1)/60)

@printf("\n")
