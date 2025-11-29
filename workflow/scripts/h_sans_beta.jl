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

h_norms1 = []
h_norms2 = []
h_norms3 = []

hrange =  10.0 .^ (-5:1:4)
hnumber = length(hrange)
h_norms = Array{Float64}(undef, hnumber, 3)
for (i,h) in enumerate(hrange)
    gk_values = generate_gk_values(function_names, calc_dyadics_over_beta, h)

    derivatives = generate_numderiv_dyadic(gk_values)

    Sans_total = Sans_beta + Sans_beta2 + Sans_beta3

    Sans_norm = norm(Sans_total, 2)

    h_diff1 = Sans_beta - derivatives[1]
    h_diff2 = Sans_beta2 - derivatives[2]
    h_diff3 = Sans_beta3 - derivatives[3]

    h_norm1 = norm(h_diff1, 2)/Sans_norm
    h_norm2 = norm(h_diff2, 2)/Sans_norm
    h_norm3 = norm(h_diff3, 2)/Sans_norm
    h_norms[i,1] = h_norm1
    h_norms[i,2] = h_norm2
    h_norms[i,3] = h_norm3
end

writedlm("$(output_dir)/h_norms.txt", h_norms)