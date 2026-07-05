# Emit a reference evaluation case for the Fortran gg field/potential kernel.
using GeneralizedGradients
using Printf, Random

const OUT = ARGS[1]

akeys = Set{Tuple{Int,Int}}(); bkeys = Set{Tuple{Int,Int}}(); bskeys = Set{Int}()
for cn in ("Bx","By","Bs","Ax","Ay","As")
  for (t,S) in ((getfield(GeneralizedGradients,Symbol("$(cn)_a")),akeys),
                (getfield(GeneralizedGradients,Symbol("$(cn)_b")),bkeys))
    for k in keys(t); push!(S,k); end
  end
  for k in keys(getfield(GeneralizedGradients,Symbol("$(cn)_bs"))); push!(bskeys,k); end
end

Random.seed!(2024)
aval  = Dict(k => randn() for k in akeys)
bval  = Dict(k => randn() for k in bkeys)
bsval = Dict(k => randn() for k in bskeys)

fit = GGCoefs()
fit.z_base = [0.0]; fit.g_ref = 0.37; fit.origin = [0.0,0.0]; fit.dz_grid = 0.1
fit.m_max = maximum(m for (_,m) in akeys)
fit.a  = Dict(k=>[v] for (k,v) in aval)
fit.b  = Dict(k=>[v] for (k,v) in bval)
fit.bs = Dict(k=>[v] for (k,v) in bsval)

x, y = 0.011, -0.007
B, A, dA = GeneralizedGradients.field_and_potential_evaluate(fit, 1, x, y)

# gval list: kind(1=a,2=b,3=bs) n m value. Emit up to m+1 (bumped orders present as
# separate keys already, value 0 if absent).
entries = Tuple{Int,Int,Int,Float64}[]
for ((n,m),v) in aval;  push!(entries,(1,n,m,v)); end
for ((n,m),v) in bval;  push!(entries,(2,n,m,v)); end
for (m,v) in bsval;     push!(entries,(3,0,m,v)); end

open(OUT,"w") do io
  @printf(io, "%.17e %.17e %.17e\n", fit.g_ref, x, y)
  println(io, length(entries))
  for (kind,n,m,v) in entries
    @printf(io, "%d %d %d %.17e\n", kind, n, m, v)
  end
  @printf(io, "%.17e %.17e %.17e\n", B[1], B[2], B[3])
  @printf(io, "%.17e %.17e %.17e\n", A[1], A[2], A[3])
  for i in 1:3
    @printf(io, "%.17e %.17e %.17e\n", dA[i,1], dA[i,2], dA[i,3])
  end
end
println("wrote $OUT : ", length(entries), " gval entries")
