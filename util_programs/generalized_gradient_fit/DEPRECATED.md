# DEPRECATED

This program (`generalized_gradient_fit`) fits the **obsolete** straight-coordinate
generalized-gradient `gen_grad_map` format and does not support a curved
reference frame (issue #2086).

It is superseded by the `gg_fit` function of the
[GeneralizedGradients.jl](https://github.com/bmad-sim/GeneralizedGradients.jl)
package, which fits the current curved-coordinate `gen_gradients` formalism used
by Bmad. Use that instead.

Existing lattice files in the old `gen_grad_map` format can be converted to the
new `gen_gradients` format with
`util_programs/gen_gradients/gen_grad_map_to_gen_gradients.py`.

This directory is retained for reference only and is no longer maintained.
