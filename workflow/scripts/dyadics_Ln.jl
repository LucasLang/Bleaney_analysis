using ParaMag
using LinearAlgebra
using Printf
using DelimitedFiles
using Serialization

include("functions.jl")


# Variables and arguments
Tmin = parse(Float64, ARGS[1]) 
Tmax = parse(Float64, ARGS[2]) 
Tinterval = parse(Float64, ARGS[3]) 
h = readdlm(ARGS[4])[1]
finitediff_functions_file = ARGS[5]
include("../"*finitediff_functions_file)  # file is provided as relative path; have to move up one directory from "scripts"
Ln = ARGS[6]
output_file = ARGS[7]

# Get the directory of the script and construct the file path
script_dir = @__DIR__
data_dir = joinpath(script_dir, "..", "data")
param_file = joinpath(data_dir, "Bkq_$(Ln)_real")
      

# Defining shparam
shparam = ParaMag.SHParam_lanthanoid(param_file, Ln)

all_diff_norms_matrix = run_dyadics(shparam)
writedlm(output_file, all_diff_norms_matrix)
