using GeneralizedGradients
# Single-plane fit with the exact GG values used in the smoke-test lattice.
fit = GGCoefs()
fit.z_base = [0.0]; fit.g_ref = 0.37; fit.origin = [0.0,0.0]; fit.dz_grid = 0.1; fit.m_max = 2
fit.a  = Dict((2,0)=>[0.3], (2,1)=>[-0.2], (2,2)=>[0.1])
fit.b  = Dict((1,0)=>[0.5], (1,1)=>[0.1], (1,2)=>[-0.2])
fit.bs = Dict{Int,Vector{Float64}}()
x, y = 0.011, -0.007
B, A, dA = GeneralizedGradients.field_and_potential_evaluate(fit, 1, x, y)
using Printf
@printf("B = %.14e %.14e %.14e\n", B[1], B[2], B[3])
@printf("A = %.14e %.14e %.14e\n", A[1], A[2], A[3])
