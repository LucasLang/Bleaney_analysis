using Printf
using LinearAlgebra

# Modified dyadics functions from ParaMag
function calc_dyadics_beta(sh, β)
    return calc_dyadic(sh, 1 / (ParaMag.kB * β))
end

function calc_dyadics_over_beta(sh, β)
    return calc_dyadics_beta(sh, β) * 1/β
end

function generate_gk_values(sh, function_names, h::Float64)
    gk_values = []  # To store the function results and preserve order

    # Loop over the function names
    for fn_name in function_names
        approximation_result = fn_name(β -> calc_dyadics_over_beta(sh, β), h)

        # Store the result
	push!(gk_values, approximation_result)
    end
    
    return gk_values
end

# this only gives the derivative (temperature-independent)
function generate_numderiv_dyadic(gk_values)
    n_order = length(gk_values)  # Determine order based on number of elements in gk_values

    # Generate derivative_dyadic_k for each k starting from 1 (derivative_dyadic_1 = g1, derivative_dyadic_2 = 2 * g2, etc.)
    derivatives = [k * gk_values[k] for k in 1:n_order]

    return derivatives
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

function run_dyadics(shparam)
    sh = ParaMag.SpinHamiltonian(shparam)

    gk_values = generate_gk_values(sh, function_names, h)

    all_diff_norms = []
    for T in Tmin:Tinterval:Tmax
        # define beta
        β = 1 / (ParaMag.kB * T)

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
        push!(all_diff_norms, diff_norms)
    end
    Tnumber = length(all_diff_norms)
    all_diff_norms_matrix = [all_diff_norms[Tindex][order] for Tindex in 1:Tnumber, order in 1:6]

    writedlm(joinpath(output_dir, "diff_norms.txt"), all_diff_norms_matrix)

end
