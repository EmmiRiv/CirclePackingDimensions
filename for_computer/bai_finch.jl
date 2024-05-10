### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 236bdd31-4496-48e0-be58-7319738fc1b6
using LinearAlgebra

using Printf

# ╔═╡ 9144b032-ed27-4910-bf0e-b26f1bf3b8ff
r_mat = (1/12)*[[-6+8im,11im],[4im,-6-8im]]

# ╔═╡ 6a4be131-7112-4a03-b30a-e6d1fe43e0f1
r_mat2 = [[-(7 + 24im)/36, -121/144], [-1/9, (-7 + 24im)/36]]

# ╔═╡ 1041df01-3192-4f28-a0d9-77733cbb6139
g_mat = (1/12)*[[2-2im,-1-5im],[4-4im,-2-10im]]

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

# compute oftype(x, y)^p efficiently, choosing the correct branch cut

function pow_oftype(x::Complex, y::Real, p::Real)
    if p >= 0
        # note: this will never be called for y < 0,
        # which would throw an error for non-integer p here
        return oftype(x, y^p)
    else
        yp = y^-p # use real power for efficiency
        return oftype(x, Complex(yp, -zero(yp))) # get correct sign of zero!
    end
end
pow_oftype(x, y, p) = oftype(x, y)^p
pow_oftype(x::Complex, y::Real, p::Complex) = oftype(x, y^p)
# Helper macro for polygamma(m, z):
#   Evaluate p[1]*c[1] + x*p[2]*c[2] + x^2*p[3]*c[3] + ...
#   where c[1] = m + 1
#         c[k] = c[k-1] * (2k+m-1)*(2k+m-2) / ((2k-1)*(2k-2)) = c[k-1] * d[k]
#         i.e. d[k] = c[k]/c[k-1] = (2k+m-1)*(2k+m-2) / ((2k-1)*(2k-2))
#   by a modified version of Horner's rule:
#      c[1] * (p[1] + d[2]*x * (p[2] + d[3]*x * (p[3] + ...))).
# The entries of p must be literal constants and there must be > 1 of them.
macro pg_horner(x, m, p...)
    k = length(p)
    me = esc(m)
    xe = esc(x)
    ex = :(($me + $(2k-1)) * ($me + $(2k-2)) * $(p[end]/((2k-1)*(2k-2))))
    for k = length(p)-1:-1:2
        cdiv = 1 / ((2k-1)*(2k-2))
        ex = :(($cdiv * ($me + $(2k-1)) * ($me + $(2k-2))) *
               ($(p[k]) + $xe * $ex))
    end
    :(($me + 1) * ($(p[1]) + $xe * $ex))
end

"""
    zeta(s, z)

Generalized zeta function defined by
```math
\\zeta(s, z)=\\sum_{k=0}^\\infty \\frac{1}{((k+z)^2)^{s/2}},
```
where any term with ``k+z=0`` is excluded.  For ``\\Re z > 0``,
this definition is equivalent to the Hurwitz zeta function
``\\sum_{k=0}^\\infty (k+z)^{-s}``.

The Riemann zeta function is recovered as ``\\zeta(s)=\\zeta(s,1)``.

External links: [Riemann zeta function](https://en.wikipedia.org/wiki/Riemann_zeta_function), [Hurwitz zeta function](https://en.wikipedia.org/wiki/Hurwitz_zeta_function)
"""

function zeta(s,z)

    # handle NaN cases
    if isnan(s) || isnan(z)
        return T <: Real ? NaN : NaN + NaN*im
    end

    x = real(z)

    # annoying s = Inf case:
    if !isfinite(s)
        if real(s) == Inf
            if x > 1 || (x >= 0.5 ? abs(z) > 1 : abs(z - round(x)) > 1)
                return zero(T) # distance to poles is > 1
            end
            x > 0 && isreal(z) && isreal(s) && return T(Inf)
        end
        throw(DomainError(s, "`s` must be finite."))  # nothing clever to return
    end

    m = s - 1
    h = 0.0

    # Algorithm is just the m-th derivative of digamma formula above,
    # with a modified cutoff of the final asymptotic expansion.

    # Note: we multiply by -(-1)^m m! in polygamma below, so this factor is
    #       pulled out of all of our derivatives.

    cutoff = 7 + real(m) + abs(imag(m)) # TODO: this cutoff is too conservative?
    if x < cutoff
        # shift using recurrence formula
        xf = floor(x)
        nx = Int(xf)
        n = ceil(Int, cutoff - nx)
        minus_s = -s
        if nx < 0 # x < 0
            # need to use (-z)^(-s) recurrence to be correct for real z < 0
            # [the general form of the recurrence term is (z^2)^(-s/2)]
            minus_z = -z
            h += pow_oftype(h, minus_z, minus_s) # v = 0 term
            if xf != z
                h += pow_oftype(h, z - nx, minus_s)
                # real(z - nx) > 0, so use correct branch cut
                # otherwise, if xf==z, then the definition skips this term
            end
            # do loop in different order, depending on the sign of s,
            # so that we are looping from largest to smallest summands and
            # can halt the loop early if possible; see issue #15946
            # FIXME: still slow for small m, large Im(z)
            if real(s) > 0
                for v in -nx-1:-1:1
                    h0= h
                    h += pow_oftype(h, minus_z - v, minus_s)
                    h == h0 && break # prevent long loop for large -x > 0
                end
            else
                for v in 1:-nx-1
                    h0= h
                    h += pow_oftype(h, minus_z - v, minus_s)
                    h == h0 && break # prevent long loop for large -x > 0
                end
            end
        else # x ≥ 0 && z != 0
            h += pow_oftype(h, z, minus_s)
        end
        # loop order depends on sign of s, as above
        if real(s) > 0
            for v in max(1,1-nx):n-1
                h0= h
                h += pow_oftype(h, z + v, minus_s)
                h == h0 && break # prevent long loop for large m
            end
        else
            for v in n-1:-1:max(1,1-nx)
                h0= h
                h += pow_oftype(h, z + v, minus_s)
                h == h0 && break # prevent long loop for large m
            end
        end
        z += n
    end

    t = inv(z)
    w = isa(t, Real) ? conj(oftype(h, t))^m : oftype(h, t)^m
    h += w * (inv(m) + 0.5*t)

    t *= t # 1/z^2
    h += w*t * @pg_horner(t,m,0.08333333333333333,-0.008333333333333333,0.003968253968253968,-0.004166666666666667,0.007575757575757576,-0.021092796092796094,0.08333333333333333,-0.4432598039215686,3.0539543302701198)

    return h
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
    matr = 2*real(gen_approx_k0(q,1) + gen_approx_F(q,2))
    return hcat(matr...)
end

# ╔═╡ 456f05a9-f100-4804-8350-483dc9e249e4
function power_method(q)
	matr = gen_Lq(q)
    vec = ones(BigFloat, (Nc+1)^2)
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
    @printf("power method iterations: %d\n", count)
    @printf("simulated circles: %d\n",(2*k0)^count)
    return current_val
end

function secant(x0,y0,x1,y1,z)
	return x0 - (y0-z) * ((x1-x0)/(y1-y0))
end

# ╔═╡ a0b1115e-d614-4595-84fe-8fb4de91ec39
function secant_method(z,x1,x2,e,its)
	k1 = x1
	k2 = x2
	y1 = power_method(k1)
	y2 = power_method(k2)
	count = 1
	while abs(y1-z)>e && count<its
		k3 = secant(k1,y1,k2,y2,z)
		k1 = k2
		y1 = y2
		k2 = k3
		y2 = power_method(k2)
		count += 1
	end
	return count,k1,y1
end

t1 = time()
k0 = parse(Float64, ARGS[1])

# ╔═╡ aa71d4ec-94aa-4b43-a6cb-1361948d9dc7
Nc = 20

# ╔═╡ 6a94dbda-1ee4-484c-a7cf-2c58aba99ef5
Lc = 20

@printf("k0: %d\n", k0)

n,q,v = secant_method(1,1.3,1.31,1e-6,100)

@printf("dimension: %f, eigenval: %f, iterations: %d\n", q, v, n)


@printf("time: %f minutes\n", (time()-t1)/60)

@printf("\n")
