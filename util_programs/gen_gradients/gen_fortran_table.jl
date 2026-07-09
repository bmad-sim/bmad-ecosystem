# Generate a Bmad Fortran module holding the gg coefficient table
# (field + vector potential, with explicit g_ref^k dependence), transcribed
# from GeneralizedGradients.jl tables/gg_coef_table.jl.
#
# Usage: julia gen_fortran_table.jl <gg_coef_table.jl> <out.f90>

using Printf

if length(ARGS) < 2
  error("usage: julia gen_fortran_table.jl <path/to/gg_coef_table.jl> <out/gg_coef_table_mod.f90>")
end
const TABLE = ARGS[1]
const OUT   = ARGS[2]

include(TABLE)

# Component index: 1=Bx 2=By 3=Bs 4=Ax 5=Ay 6=As
# Kind index:      1=a  2=b  3=bs
const COMPS = ["Bx","By","Bs","Ax","Ay","As"]
const KINDS = Dict("a"=>1, "b"=>2, "bs"=>3)

# Collect all terms into a flat list of (comp,kind,n,m,p,q,k,coef).
terms = Tuple{Int,Int,Int,Int,Int,Int,Int,Float64}[]

for (ci, cname) in enumerate(COMPS)
  # a and b tables are keyed by (n,m); bs by m (n treated as 0).
  for kname in ("a","b")
    d = getfield(Main, Symbol("$(cname)_$(kname)"))
    for ((n,m), tlist) in d
      for (c,p,q,k) in tlist
        push!(terms, (ci, KINDS[kname], n, m, p, q, k, float(c)))
      end
    end
  end
  dbs = getfield(Main, Symbol("$(cname)_bs"))
  for (m, tlist) in dbs
    for (c,p,q,k) in tlist
      push!(terms, (ci, 3, 0, m, p, q, k, float(c)))
    end
  end
end

# Sort for reproducible output.
sort!(terms)
nterm = length(terms)

# Sizing constants.
maxn = maximum(t[3] for t in terms)
maxm = maximum(t[4] for t in terms)
maxpq = maximum(max(t[5],t[6]) for t in terms)
maxk = maximum(t[7] for t in terms)

open(OUT, "w") do io
  println(io, "!+")
  println(io, "! Module gg_coef_table_mod")
  println(io, "!")
  println(io, "! GENERATED FILE -- do not edit by hand.")
  println(io, "! Produced from GeneralizedGradients.jl tables/gg_coef_table.jl by")
  println(io, "! util_programs/gen_gradients/gen_fortran_table.jl .")
  println(io, "!")
  println(io, "! Holds the generalized-gradient coefficient table used to evaluate the")
  println(io, "! magnetic field B = (Bx,By,Bs) and vector potential A = (Ax,Ay,As) in a")
  println(io, "! straight OR constant-curvature (g_ref = 1/rho) reference frame.")
  println(io, "!")
  println(io, "! For output component `comp` (1=Bx 2=By 3=Bs 4=Ax 5=Ay 6=As) and gg source")
  println(io, "! `kind` (1=a skew, 2=b normal, 3=bs solenoid), each table entry contributes")
  println(io, "!   comp += coef * g_ref**k * x**p * y**q * G(kind,n,m)")
  println(io, "! where G(a,n,m)=d^m a_n/ds^m, G(b,n,m)=d^m b_n/ds^m, G(bs,0,m)=d^m b_s/ds^m.")
  println(io, "!-")
  println(io, "")
  println(io, "module gg_coef_table_mod")
  println(io, "")
  println(io, "use precision_def")
  println(io, "")
  println(io, "implicit none")
  println(io, "")
  println(io, "integer, parameter :: gg_n_coef_term\$ = $nterm")
  println(io, "integer, parameter :: gg_coef_max_n\$ = $maxn, gg_coef_max_m\$ = $maxm")
  println(io, "integer, parameter :: gg_coef_max_pq\$ = $maxpq, gg_coef_max_k\$ = $maxk")
  println(io, "")
  println(io, "type gg_coef_term_struct")
  println(io, "  integer :: comp = 0    ! 1=Bx 2=By 3=Bs 4=Ax 5=Ay 6=As")
  println(io, "  integer :: kind = 0    ! 1=a 2=b 3=bs")
  println(io, "  integer :: n = 0       ! Harmonic index (0 for bs)")
  println(io, "  integer :: m = 0       ! s-derivative order")
  println(io, "  integer :: p = 0       ! x power")
  println(io, "  integer :: q = 0       ! y power")
  println(io, "  integer :: k = 0       ! g_ref power")
  println(io, "  real(rp) :: coef = 0")
  println(io, "end type")
  println(io, "")
  println(io, "type (gg_coef_term_struct), protected, save :: gg_coef_table(gg_n_coef_term\$)")
  println(io, "logical, private, save :: gg_table_initted = .false.")
  println(io, "")
  println(io, "contains")
  println(io, "")
  println(io, "!+")
  println(io, "! Subroutine gg_coef_table_init ()")
  println(io, "!")
  println(io, "! Populate the module-level gg_coef_table on first call (idempotent).")
  println(io, "!-")
  println(io, "")
  println(io, "subroutine gg_coef_table_init ()")
  println(io, "if (gg_table_initted) return")

  # Split the assignments into blocks to keep each subroutine small.
  blocksize = 600
  nblock = cld(nterm, blocksize)
  for b in 1:nblock
    println(io, "call gg_set_block_$(lpad(b,3,'0')) ()")
  end
  println(io, "gg_table_initted = .true.")
  println(io, "end subroutine gg_coef_table_init")
  println(io, "")

  for b in 1:nblock
    i0 = (b-1)*blocksize + 1
    i1 = min(b*blocksize, nterm)
    println(io, "subroutine gg_set_block_$(lpad(b,3,'0')) ()")
    for i in i0:i1
      (comp,kind,n,m,p,q,k,coef) = terms[i]
      cstr = coef == floor(coef) && abs(coef) < 1e15 ?
             "$(Int(coef)).0_rp" : @sprintf("%.17e_rp", coef)
      # Fortran doesn't accept 'e' exponent with _rp kind suffix cleanly for
      # some compilers; keep plain decimal via %.17e then it's fine as real.
      println(io, "gg_coef_table($i) = gg_coef_term_struct($comp,$kind,$n,$m,$p,$q,$k, $cstr)")
    end
    println(io, "end subroutine gg_set_block_$(lpad(b,3,'0'))")
    println(io, "")
  end

  println(io, "end module gg_coef_table_mod")
end

println("Wrote $OUT : $nterm terms, maxn=$maxn maxm=$maxm maxpq=$maxpq maxk=$maxk, $(cld(nterm,600)) blocks")
