using Printf
using LinearAlgebra

# Modified dyadics functions from ParaMag
function calc_dyadics_beta(sh, β)
    return calc_dyadic(sh, 1 / (ParaMag.kB * β))
end

function generate_derivatives(sh, function_names, h::Float64)
    gk_values = []  # To store the function results and preserve order

    # Loop over the function names
    for fn_name in function_names
        approximation_result = fn_name(β -> calc_dyadics_beta(sh, β), h)

        # Store the result
	push!(gk_values, approximation_result)
    end
    
    return gk_values[2:end]   # exclude the first value, which is just the dyadic at beta=0 (which is zero)
end

#########################################################################################
# Output settings

script_dir = @__DIR__
output_dir = joinpath(script_dir, "..", "results")

# Helper function to format and write matrices to a file
function write_matrix(file, name, matrix)
    println(file, "\n$name:")
    for i in 1:size(matrix, 1)
        for j in 1:size(matrix, 2)
            @printf(file, " %12.6e ", matrix[i, j])
        end
        println(file)
    end
end

function run_dyadics(shparam, temperatures)
    sh = ParaMag.SpinHamiltonian(shparam)

    num_derivatives = generate_derivatives(sh, function_names, h)

    all_diff_norms = []
    for T in temperatures
        # define beta
        β = 1 / (ParaMag.kB * T)

        # exact Dyadics
        exact_dyadic = calc_dyadic(sh, T)

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
        push!(all_diff_norms, diff_norms)
    end
    Tnumber = length(all_diff_norms)
    all_diff_norms_matrix = [all_diff_norms[Tindex][order] for Tindex in 1:Tnumber, order in 1:6]
    return all_diff_norms_matrix
end

function run_hnorms(shparam, h_minexp, h_maxexp)
    sh = ParaMag.SpinHamiltonian(shparam)

    Sans_beta = ParaMag.JJbeta(shparam)
    Sans_beta2 = ParaMag.JJbeta2(shparam)
    Sans_beta3 = ParaMag.JJbeta3(shparam)

    h_norms1 = []
    h_norms2 = []
    h_norms3 = []

    exprange = h_minexp:1:h_maxexp
    hrange = 10.0 .^ exprange
    hnumber = length(hrange)
    h_norms = Array{Float64}(undef, hnumber, 3)
    for (i,h) in enumerate(hrange)
        derivatives = generate_derivatives(sh, function_names, h)

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
    return h_norms
end
