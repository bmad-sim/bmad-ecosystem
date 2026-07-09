# gen_gradients tools

Tooling for the curved-coordinate generalized-gradient (GG) formalism
(issue #2086), which ports the field + vector-potential coefficient tables from
[GeneralizedGradients.jl](https://github.com/bmad-sim/GeneralizedGradients.jl)
into Bmad.

## `gen_fortran_table.jl`

Transcribes the GG coefficient table (`tables/gg_coef_table.jl` in
GeneralizedGradients.jl) into the Bmad Fortran module
`bmad/modules/gg_coef_table_mod.f90`.

```
julia gen_fortran_table.jl <path/to/GeneralizedGradients.jl/tables/gg_coef_table.jl> \
                           <path/to/bmad-ecosystem/bmad/modules/gg_coef_table_mod.f90>
```

The generated module holds, for every output component (Bx,By,Bs,Ax,Ay,As) and
GG source (`a` skew, `b` normal, `bs` solenoid), the list of monomial terms
`coef * g_ref**k * x**p * y**q * G(kind,n,m)` where `g_ref = 1/rho` is the
constant curvature of the reference frame. The field and vector potential are
evaluated by summing these terms into an (x,y) coefficient array and then
evaluating the polynomial (see `gg_coef_table_mod` header). `d/ds` derivatives
are obtained by bumping the GG derivative order (`a(n,m) -> a(n,m+1)`).

This only needs to be re-run if the derivative-order cutoff in
GeneralizedGradients.jl is extended (currently orders to 12, giving 7944 terms).

## Validation scripts

- `validate_table.jl <gg_coef_table_mod.f90>` — parses the generated Fortran
  module and checks every coefficient round-trips exactly against the Julia
  source dicts.
- `func_validate.jl <gg_coef_table_mod.f90>` — builds a random single-plane
  `GGCoefs` with a curved frame and checks that the table-driven eval formula
  (the one the Bmad field-calc implements) reproduces
  `GeneralizedGradients.field_and_potential_evaluate` for B, A, and dA to
  machine precision. Run with `julia --project=<GeneralizedGradients.jl> ...`.
- `emit_ref.jl <gg_ref_case.txt>` — writes a reference evaluation case (curved
  frame, random GG values, plus B/A/dA from the package).
- `test_kernel.f90` — standalone Fortran driver that reads `gg_ref_case.txt`,
  runs `gg_field_potential_calc` (bmad/modules/gen_gradients_mod.f90), and checks
  it against the reference. Build alongside `precision_def`, `gg_coef_table_mod`,
  `gen_gradients_mod` with `gfortran -ffree-line-length-none -fdollar-ok`.
