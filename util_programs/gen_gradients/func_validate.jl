# Functional validation: reproduce GeneralizedGradients.field_and_potential_evaluate
# using the PARSED Fortran table + the accum/polyval formula the Bmad field-calc
# will implement. Curved frame (g_ref != 0).
using GeneralizedGradients
using GeneralizedGradients: Bx_a, By_a, Bs_a, Bx_b, By_b, Bs_b
using Printf
using Random

const F90 = ARGS[1]

# --- parse the generated Fortran table ---
terms = Tuple{Int,Int,Int,Int,Int,Int,Int,Float64}[]
rx = r"gg_coef_term_struct\(\s*(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+),\s*([^)]+)\)"
for line in eachline(F90)
  mm = match(rx, line); mm===nothing && continue
  v = mm.captures
  push!(terms, (parse(Int,v[1]),parse(Int,v[2]),parse(Int,v[3]),parse(Int,v[4]),
                parse(Int,v[5]),parse(Int,v[6]),parse(Int,v[7]),
                parse(Float64, replace(strip(v[8]),"_rp"=>""))))
end

# --- gather (n,m) keys for a/b and m for bs across all component tables ---
akeys = Set{Tuple{Int,Int}}(); bkeys = Set{Tuple{Int,Int}}(); bskeys = Set{Int}()
for cn in ("Bx","By","Bs","Ax","Ay","As")
  for (t,S) in ((getfield(GeneralizedGradients,Symbol("$(cn)_a")),akeys),
                (getfield(GeneralizedGradients,Symbol("$(cn)_b")),bkeys))
    for k in keys(t); push!(S,k); end
  end
  for k in keys(getfield(GeneralizedGradients,Symbol("$(cn)_bs"))); push!(bskeys,k); end
end

Random.seed!(12345)
aval = Dict(k => randn() for k in akeys)
bval = Dict(k => randn() for k in bkeys)
bsval = Dict(k => randn() for k in bskeys)

# --- build a single-plane GGCoefs ---
fit = GGCoefs()
fit.z_base = [0.0]
fit.g_ref  = 0.37          # curved frame
fit.origin = [0.0, 0.0]
fit.dz_grid = 0.1
fit.m_max  = maximum(m for (_,m) in akeys)
fit.a  = Dict(k => [v] for (k,v) in aval)
fit.b  = Dict(k => [v] for (k,v) in bval)
fit.bs = Dict(k => [v] for (k,v) in bsval)

x, y = 0.011, -0.007
Bref, Aref, dAref = GeneralizedGradients.field_and_potential_evaluate(fit, 1, x, y)

# --- my table-driven formula (what Fortran will do) ---
NP = 16
function comp_arrays(g_ref, aval, bval, bsval)
  K = zeros(6, NP, NP)
  for (comp,kind,n,m,p,q,k,coef) in terms
    g = kind==1 ? get(aval,(n,m),0.0) : kind==2 ? get(bval,(n,m),0.0) : get(bsval,m,0.0)
    g == 0.0 && continue
    K[comp,p+1,q+1] += coef * (k==0 ? 1.0 : g_ref^k) * g
  end
  return K
end
polyval(K,c,x,y) = sum(K[c,i+1,j+1]*x^i*y^j for i in 0:NP-1, j in 0:NP-1)
polyval_dx(K,c,x,y) = sum(i==0 ? 0.0 : K[c,i+1,j+1]*i*x^(i-1)*y^j for i in 0:NP-1, j in 0:NP-1)
polyval_dy(K,c,x,y) = sum(j==0 ? 0.0 : K[c,i+1,j+1]*j*x^i*y^(j-1) for i in 0:NP-1, j in 0:NP-1)

K  = comp_arrays(fit.g_ref, aval, bval, bsval)
# bumped orders for d/ds:  a(n,m)->a(n,m+1)
avalp  = Dict(k => get(aval,(k[1],k[2]+1),0.0) for k in akeys)
bvalp  = Dict(k => get(bval,(k[1],k[2]+1),0.0) for k in bkeys)
bsvalp = Dict(m => get(bsval,m+1,0.0) for m in bskeys)
Kp = comp_arrays(fit.g_ref, avalp, bvalp, bsvalp)

Bmine = [polyval(K,1,x,y), polyval(K,2,x,y), polyval(K,3,x,y)]
Amine = [polyval(K,4,x,y), polyval(K,5,x,y), polyval(K,6,x,y)]
dAmine = [polyval_dx(K,4,x,y) polyval_dy(K,4,x,y) polyval(Kp,4,x,y);
          polyval_dx(K,5,x,y) polyval_dy(K,5,x,y) polyval(Kp,5,x,y);
          polyval_dx(K,6,x,y) polyval_dy(K,6,x,y) polyval(Kp,6,x,y)]

relerr(a,b) = maximum(abs.(a .- b) ./ max.(abs.(a), 1e-30))
@printf("B  max rel err = %.3e\n", relerr(Bref, Bmine))
@printf("A  max rel err = %.3e\n", relerr(Aref, Amine))
@printf("dA max rel err = %.3e\n", relerr(dAref, dAmine))
ok = relerr(Bref,Bmine)<1e-10 && relerr(Aref,Amine)<1e-10 && relerr(dAref,dAmine)<1e-10
println(ok ? "FUNCTIONAL VALIDATION PASSED" : "FUNCTIONAL VALIDATION FAILED")
