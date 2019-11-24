import math
local_funcs = {
    "acos": math.acos, 
    "acosh": math.acosh, 
    "asin": math.asin, 
    "asinh": math.asinh, 
    "atan": math.atan, 
    "atan2": math.atan2, 
    "atanh": math.atanh, 
    "ceil": math.ceil, 
    "cos": math.cos, 
    "cosh": math.cosh, 
    "degrees": math.degrees,
    "e_log": math.e,
    "erf": math.erf, 
    "erfc": math.erfc, 
    "exp": math.exp, 
    "factorial": math.factorial, 
    "floor": math.floor, 
    "gamma": math.gamma, 
    "log": math.log, 
    "log10": math.log10, 
    "pi": math.pi, 
    "pow": math.pow, 
    "radians": math.radians, 
    "sin": math.sin, 
    "sinh": math.sinh, 
    "sqrt": math.sqrt,
    "tan": math.tan,
    "tanh": math.tanh}

def eval_expr(expr):
  if "__" in expr or "=" in expr: raise  # Disallow unsafe evaluations.
  return eval(expr, local_funcs)
