using ParaMag
using LinearAlgebra
using Printf
using DelimitedFiles
using Serialization

include("functions.jl")

# Variables and arguments
Dtensor_file = ARGS[1]
S = parse(Float64, ARGS[2])
temperatures_file = ARGS[3]
h = readdlm(ARGS[4])[1]
finitediff_functions_file = ARGS[5]
output_file = ARGS[6]

include("../"*finitediff_functions_file)  # file is provided as relative path; have to move up one directory from "scripts"
gtensor = 2*Matrix(1.0I, 3, 3)

# Defining shparam
mult = Int64(2 * S + 1)
Dtensor = readdlm(Dtensor_file)
shparam = ParaMag.SHParam(mult, gtensor, Dtensor)

temperatures = readdlm(temperatures_file)
all_diff_norms_matrix = run_dyadics(shparam, temperatures)
writedlm(output_file, all_diff_norms_matrix)