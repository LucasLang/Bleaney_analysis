using ParaMag
using LinearAlgebra
using Printf
using DelimitedFiles
include("functions.jl")

h_minexp = parse(Int64, ARGS[1])
h_maxexp = parse(Int64, ARGS[2])
finitediff_functions_file = ARGS[3]
include("../"*finitediff_functions_file)  # file is provided as relative path; have to move up one directory from "scripts"
Dtensor_file = ARGS[4]
S = parse(Float64, ARGS[5])
outfile = ARGS[6]
         
gtensor = 2*Matrix(1.0I, 3, 3)
mult = Int64(2 * S + 1)
Dtensor = readdlm(Dtensor_file)
shparam = ParaMag.SHParam(mult, gtensor, Dtensor)

h_norms = run_hnorms(shparam, h_minexp, h_maxexp)
writedlm(outfile, h_norms)