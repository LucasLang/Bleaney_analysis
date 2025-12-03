using ParaMag
using LinearAlgebra
using Printf
using DelimitedFiles
using Serialization

include("functions.jl")

# Input settings
# Get the directory of the script and construct the file path
script_dir = @__DIR__
data_dir = joinpath(script_dir, "..", "data")
output_dir = joinpath(script_dir, "..", "results")
file_path = joinpath(data_dir, "spin_dyadics_input.txt")

# Read and extract Bkq file name
function extract_ln(folder)
    files = readdir(folder)

    for file in files
        m2 = match(r"Bkq_(.+?)_real", file)
        if m2 !== nothing  
            return String(m2[1]) 
        end
    end

    return nothing
end

Ln = extract_ln(data_dir)
fallback_file = joinpath(data_dir, "Bkq_$(Ln)_real")
      
# Variables and arguments
Tmin = parse(Float64, ARGS[1]) 
Tmax = parse(Float64, ARGS[2]) 
Tinterval = parse(Float64, ARGS[3]) 
h = readdlm(ARGS[4])[1]
finitediff_functions_file = ARGS[5]
include("../"*finitediff_functions_file)  # file is provided as relative path; have to move up one directory from "scripts"
gtensor = 2*Matrix(1.0I, 3, 3)

# Defining shparam
if isfile(file_path) && isfile(fallback_file)
    error("Both the primary file and the fallback file are available. Please remove one to proceed.")
elseif isfile(file_path)
    S_path = joinpath(data_dir, "S.txt")
    S_file = open(S_path, "r")
    S = parse(Float64, readline(S_file))
    mult = Int64(2 * S + 1)
    Dtensor = readdlm(file_path)
    shparam = ParaMag.SHParam(mult, gtensor, Dtensor)
elseif isfile(fallback_file)
    shparam = ParaMag.SHParam_lanthanoid(fallback_file, Ln)
else
    error("Both the primary and fallback files are missing. Cannot proceed.")
end

run_dyadics(shparam)
