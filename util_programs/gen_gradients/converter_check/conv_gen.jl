using GeneralizedGradients
using Printf
outbase = ARGS[1]

npl = 3; dz = 0.1
fit = GGCoefs()
fit.z_base = collect(0.0:dz:(npl-1)*dz)
fit.g_ref = 0.0; fit.origin = [0.0,0.0]; fit.dz_grid = dz; fit.m_max = 2

import Random; Random.seed!(7)
rv() = round.(randn(npl); digits=6)
fit.a = Dict((1,0)=>rv(),(1,1)=>rv(),(1,2)=>rv(),(2,0)=>rv(),(2,1)=>rv(),(2,2)=>rv(),(3,0)=>rv(),(3,1)=>rv(),(3,2)=>rv())
fit.b = Dict((1,0)=>rv(),(1,1)=>rv(),(1,2)=>rv(),(2,0)=>rv(),(2,1)=>rv(),(2,2)=>rv(),(3,0)=>rv(),(3,1)=>rv(),(3,2)=>rv())
fit.bs = Dict(0=>rv(), 1=>rv(), 2=>rv())

# Write the old gen_grad_map block (into <outbase>_gg.bmad).
GeneralizedGradients.write_bmad_gg_fit(fit; output_base=outbase)

# Dump the original a/b/bs at plane ip=2 (z=0.1) for the recovery check.
ip = 2
open(outbase * "_orig.txt", "w") do io
  for (nm,v) in sort(collect(fit.a)); @printf(io, "a %d %d %.14e\n", nm[1], nm[2], v[ip]); end
  for (nm,v) in sort(collect(fit.b)); @printf(io, "b %d %d %.14e\n", nm[1], nm[2], v[ip]); end
  for (m,v) in sort(collect(fit.bs)); @printf(io, "bs 0 %d %.14e\n", m, v[ip]); end
end

# Reference field at that plane.
B, A, dA = GeneralizedGradients.field_and_potential_evaluate(fit, ip, 0.011, -0.007)
@printf("REF B = %.14e %.14e %.14e\n", B[1], B[2], B[3])
