# Validate the generated Fortran table against the Julia source dicts.
using Printf
# usage: julia validate_table.jl <gg_coef_table_mod.f90> <gg_coef_table.jl>
const F90   = ARGS[1]
const TABLE = ARGS[2]
include(TABLE)

const COMPS = ["Bx","By","Bs","Ax","Ay","As"]
const KINDS = Dict("a"=>1, "b"=>2, "bs"=>3)

# Rebuild the reference flat multiset exactly as the generator did.
ref = Tuple{Int,Int,Int,Int,Int,Int,Int,Float64}[]
for (ci,cname) in enumerate(COMPS)
  for kname in ("a","b")
    d = getfield(Main, Symbol("$(cname)_$(kname)"))
    for ((n,m),tl) in d, (c,p,q,k) in tl
      push!(ref, (ci,KINDS[kname],n,m,p,q,k,float(c)))
    end
  end
  for (m,tl) in getfield(Main, Symbol("$(cname)_bs")), (c,p,q,k) in tl
    push!(ref, (ci,3,0,m,p,q,k,float(c)))
  end
end
sort!(ref)

# Parse the Fortran file.
got = Tuple{Int,Int,Int,Int,Int,Int,Int,Float64}[]
rx = r"gg_coef_term_struct\(\s*(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+),(-?\d+),\s*([^)]+)\)"
for line in eachline(F90)
  m = match(rx, line)
  m === nothing && continue
  v = m.captures
  coefstr = replace(strip(v[8]), "_rp"=>"")
  push!(got, (parse(Int,v[1]),parse(Int,v[2]),parse(Int,v[3]),parse(Int,v[4]),
              parse(Int,v[5]),parse(Int,v[6]),parse(Int,v[7]),parse(Float64,coefstr)))
end
sort!(got)

@printf("ref terms = %d,  parsed terms = %d\n", length(ref), length(got))
if length(ref) != length(got)
  println("LENGTH MISMATCH"); exit(1)
end

maxrel = 0.0; nbad = 0
for (a,b) in zip(ref,got)
  if a[1:7] != b[1:7]
    println("KEY MISMATCH: ", a, " vs ", b); global nbad += 1; nbad>5 && break
    continue
  end
  ca, cb = a[8], b[8]
  rel = ca == 0 ? abs(cb) : abs(cb-ca)/abs(ca)
  global maxrel = max(maxrel, rel)
  if rel > 1e-12; println("COEF MISMATCH ", a[1:7], " ref=", ca, " got=", cb); global nbad+=1; nbad>10 && break; end
end
@printf("max relative coef error = %.3e,  nbad = %d\n", maxrel, nbad)
println(nbad==0 && maxrel < 1e-12 ? "VALIDATION PASSED" : "VALIDATION FAILED")
