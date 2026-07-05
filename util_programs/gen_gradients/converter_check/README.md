# Converter round-trip check

Validates `gen_grad_map_to_gen_gradients.py` (old azimuthal-harmonic
`gen_grad_map` -> new midplane `gen_gradients`).

1. `conv_gen.jl <base>` — build a straight (`g_ref=0`) `GGCoefs`, write the old
   `gen_grad_map` block via `write_bmad_gg_fit` (`<base>_gg.bmad`), dump the
   original `a/b/bs` at a plane (`<base>_orig.txt`), and print the reference
   field. Run with `julia --project=<GeneralizedGradients.jl>`.
2. `gen_grad_map_to_gen_gradients.py <base>_gg.bmad -o new_block.bmad` — convert.
3. `conv_compare.py new_block.bmad <base>_orig.txt` — the recovered `a/b/bs`
   match the originals to 0.0 relative error (the inverse recursions are exact).
4. Wrapping `new_block.bmad` in an `em_field` element and evaluating with Bmad
   reproduces the `GeneralizedGradients` reference field to 14 digits.
