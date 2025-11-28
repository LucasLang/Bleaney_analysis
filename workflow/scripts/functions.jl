using Printf
using LinearAlgebra

# Modified dyadics functions from ParaMag
function calc_dyadics_beta(β)
    return calc_dyadic(sh, 1 / (ParaMag.kB * β))
end

function calc_dyadics_over_beta(β)
    return calc_dyadics_beta(β) * 1/β
end

function generate_gk_values(function_names, calc_dyadics_over_beta::Function, h::Float64)
    gk_values = []  # To store the function results and preserve order

    # Loop over the function names
    for fn_name in function_names
        approximation_result = fn_name(calc_dyadics_over_beta, h)

        # Store the result
	push!(gk_values, approximation_result)
    end
    
    return gk_values
end

# this gives a complete term in the Taylor series expansion (temperature-dependent)
function generate_dyadic_terms(β, gk_values)
    n_order = length(gk_values)  # Determine order based on number of elements in gk_values

    # Generate derivative_dyadic_k for each k starting from 1 (derivative_dyadic_1 = g1, derivative_dyadic_2 = 2 * g2, etc.)
    derivatives = [k * gk_values[k] for k in 1:n_order]

    # Generate full_beta_derivative_k for each k with factorial term and β power
    full_beta_derivatives = [
        derivatives[k] * β^k / factorial(k) for k in 1:n_order
    ]
    
    return full_beta_derivatives
end

# this only gives the derivative (temperature-independent)
function generate_dyadic_terms2(gk_values)
    n_order = length(gk_values)  # Determine order based on number of elements in gk_values

    # Generate derivative_dyadic_k for each k starting from 1 (derivative_dyadic_1 = g1, derivative_dyadic_2 = 2 * g2, etc.)
    derivatives = [k * gk_values[k] for k in 1:n_order]

    return derivatives
end

# Test
#full_beta_derivatives = generate_dyadic_terms(β, gk_values)
#println(full_beta_derivatives)

#########################################################################################
# Output settings

script_dir = @__DIR__
output_dir = joinpath(script_dir, "..", "results")

# Printing plot points into files
function write_norms_to_file()
    open(joinpath(output_dir, "diff_norm1.txt"), "a") do f1
        println(f1, diff_norm1)
    end

    open(joinpath(output_dir, "diff_norm2.txt"), "a") do f2
        println(f2, diff_norm2)
    end

    open(joinpath(output_dir, "diff_norm3.txt"), "a") do f3
        println(f3, diff_norm3)
    end

    open(joinpath(output_dir, "diff_norm4.txt"), "a") do f4
        println(f4, diff_norm4)
    end

    open(joinpath(output_dir, "diff_norm5.txt"), "a") do f5
        println(f5, diff_norm5)
    end

    open(joinpath(output_dir, "diff_norm6.txt"), "a") do f6
        println(f6, diff_norm6)
    end
end

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
