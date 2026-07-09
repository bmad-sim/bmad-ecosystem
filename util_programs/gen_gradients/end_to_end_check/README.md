# End-to-end check

Confirms the full Bmad path for the curved generalized-gradient formalism:
`bmad_parser` reads the new `gen_gradients` lattice block, `em_field_calc`
builds the GG value array and calls `gg_field_potential_calc`, and the result
matches `GeneralizedGradients.jl` exactly.

- `gg_smoke.bmad` — an `em_field` element with a curved (`g_ref = 0.37`)
  `gen_gradients` map (a `b`, `n=1` curve and an `a`, `n=2` curve).
- `gg_smoke.f90` — parses the lattice and prints `field%B` / `field%A` at a base
  plane, `(x,y) = (0.011, -0.007)`.
- `smoke_ref.jl` — the same GG values evaluated with
  `GeneralizedGradients.field_and_potential_evaluate`.

Result (both, to 14 digits):
```
B =  3.29728421958360E-03  5.02113405731874E-01 -7.04313578867046E-04
A = -2.51299243861083E-06  0.00000000000000E+00 -5.51194652714348E-03
```
