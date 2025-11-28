using ParaMag
using LinearAlgebra
using Printf
using DelimitedFiles
using Serialization

include("../results/finitediff_functions.jl")    #XXXLucasXXX: need to remove this hard-coded

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
T = parse(Float64, ARGS[1]) 
h = parse(Float64, ARGS[2])  
n_order = parse(Int, ARGS[3])
n_stencil = parse(Int, ARGS[4])
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

sh = ParaMag.SpinHamiltonian(shparam)

# Define the cache file
cache_file = joinpath(output_dir, "derivatives_results.jld2")
marker_file = joinpath(output_dir, "derivatives_done.marker")

if !isfile(marker_file)
    println("Running numerical derivatives calculator...")

    gk_values = generate_gk_values(function_names, calc_dyadics_over_beta, h)
    
    # Save precomputed results to a file
    println("Saving precomputed results to file...")
    #@save cache_file gk_values
    open(cache_file, "w") do io
        serialize(io, gk_values)
    end

    # Create the marker file
    write(marker_file, "done")
    println("Derivatives calculations complete.")
else
    println("Loading cached results...")
    #@load cache_file gk_values
    gk_values = open(cache_file, "r") do io
        deserialize(io)
    end
end

# define beta
β = 1 / (ParaMag.kB * T)
#println("β sucessfully defined!")

# exact Dyadics
exact_dyadics = calc_dyadic(sh, T)

# JJ_beta = analytical derivatives
first_deriv = ParaMag.JJbeta(shparam)
second_deriv = ParaMag.JJbeta2(shparam)
third_deriv = ParaMag.JJbeta3(shparam)

SS_beta = first_deriv * β
SS_beta2 = 1/2 * second_deriv * β^2
SS_beta3 = 1/6 * third_deriv * β^3

SS_Curie = SS_beta
SS_Bleaney = SS_Curie + SS_beta2
SS_Lucas = SS_Bleaney + SS_beta3

# Numerical derivatives calculator section
full_beta_derivatives = generate_dyadic_terms(β, gk_values)

# Plot points generation
# beta plots (analytical-numerical)

diff_beta1 = exact_dyadics - SS_Curie
diff_beta2 = exact_dyadics - SS_Bleaney
diff_beta3 = exact_dyadics - SS_Lucas
diff_beta4 = diff_beta3 - full_beta_derivatives[4]
diff_beta5 = diff_beta4 - full_beta_derivatives[5]
diff_beta6 = diff_beta5 - full_beta_derivatives[6]

exact_norm = norm(exact_dyadics, 2)

# norms
diff_norm1 = norm(diff_beta1, 2)/exact_norm
diff_norm2 = norm(diff_beta2, 2)/exact_norm
diff_norm3 = norm(diff_beta3, 2)/exact_norm
diff_norm4 = norm(diff_beta4, 2)/exact_norm
diff_norm5 = norm(diff_beta5, 2)/exact_norm
diff_norm6 = norm(diff_beta6, 2)/exact_norm

write_norms_to_file()