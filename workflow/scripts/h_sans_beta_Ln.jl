using ParaMag
using LinearAlgebra
using Printf
using DelimitedFiles
include("functions.jl")

h_minexp = parse(Int64, ARGS[1])
h_maxexp = parse(Int64, ARGS[2])
finitediff_functions_file = ARGS[3]
Ln = ARGS[4]
outfile = ARGS[5]
include("../"*finitediff_functions_file)  # file is provided as relative path; have to move up one directory from "scripts"

# Input settings
script_dir = @__DIR__
data_dir = joinpath(script_dir, "..", "data")
param_file = joinpath(data_dir, "Bkq_$(Ln)_real")

shparam = ParaMag.SHParam_lanthanoid(param_file, Ln)

h_norms = run_hnorms(shparam, h_minexp, h_maxexp)
writedlm(outfile, h_norms)