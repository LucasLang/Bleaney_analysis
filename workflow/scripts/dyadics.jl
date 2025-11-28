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
exact_dyadic = calc_dyadic(sh, T)

num_derivatives = generate_numderiv_dyadic(gk_values)
# for the first three, we use analytical derivatives, for the next three numerical
derivatives = [ParaMag.JJbeta(shparam), ParaMag.JJbeta2(shparam), ParaMag.JJbeta3(shparam),
               num_derivatives[4], num_derivatives[5], num_derivatives[6]]

Taylor_terms = [derivatives[k] * β^k / factorial(k) for k in 1:6]

# Plot points generation
# beta plots (analytical-numerical)
diff_betas = [exact_dyadic - Taylor_terms[1]]   # difference when we approximate the dyadic with just the linear (Curie) term
for order in 2:6
    push!(diff_betas, diff_betas[order-1] - Taylor_terms[order])
end

exact_norm = norm(exact_dyadic, 2)

# norms
diff_norms = [norm(diff_beta, 2)/exact_norm for diff_beta in diff_betas]

write_norms_to_file()