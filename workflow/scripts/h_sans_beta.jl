using ParaMag
using LinearAlgebra
using Printf
using DelimitedFiles
#using DataStructures
include("functions.jl")
include("../results/finitediff_functions.jl")    #XXXLucasXXX: need to remove this hard-coded

# Input settings
# Get the directory of the script and construct the file path
script_dir = @__DIR__
data_dir = joinpath(script_dir, "..", "data")
output_dir = joinpath(script_dir, "..", "results")
file_path = joinpath(data_dir, "spin_dyadics_input.txt")


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

# Constants
         
# Variables as arguments
T = parse(Float64, ARGS[1])  
h = parse(Float64, ARGS[2])  
n_order = parse(Int, ARGS[3])
n_stencil = parse(Int, ARGS[4])

gtensor = 2*Matrix(1.0I, 3, 3)

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

Sans_beta = ParaMag.JJbeta(shparam)
Sans_beta2 = ParaMag.JJbeta2(shparam)
Sans_beta3 = ParaMag.JJbeta3(shparam)

gk_values = generate_gk_values(function_names, calc_dyadics_over_beta, h)

derivatives = generate_dyadic_terms2(gk_values)

Sans_total = Sans_beta + Sans_beta2 + Sans_beta3

Sans_norm = norm(Sans_total, 2)

h_diff1 = Sans_beta - derivatives[1]
h_diff2 = Sans_beta2 - derivatives[2]
h_diff3 = Sans_beta3 - derivatives[3]

h_norm1 = norm(h_diff1, 2)/Sans_norm
h_norm2 = norm(h_diff2, 2)/Sans_norm
h_norm3 = norm(h_diff3, 2)/Sans_norm

# Printing plot points into files
function write_norms_to_file()
    open(joinpath(output_dir, "h_norm1.txt"), "a") do f1
        println(f1, h_norm1)
    end

    open(joinpath(output_dir, "h_norm2.txt"), "a") do f2
        println(f2, h_norm2)
    end

    open(joinpath(output_dir, "h_norm3.txt"), "a") do f3
        println(f3, h_norm3)
    end

    #open("h_norm4.txt", "a") do f4
        #println(f4, h_norm4)
    #end

    #open("h_norm5.txt", "a") do f5
        #println(f5, h_norm5)
    #end

    #open("h_norm6.txt", "a") do f6
        #println(f6, h_norm6)
    #end

end

write_norms_to_file()
