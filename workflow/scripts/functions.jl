using Printf
using LinearAlgebra

# Modified dyadics functions from ParaMag
function calc_dyadics_beta(β)
    return calc_dyadic(sh, 1 / (ParaMag.kB * β))
end

function calc_dyadics_over_beta(β)
    return calc_dyadics_beta(β) * 1/β
end

# Creation of Taylor's expansion and their respective coefficients
function generate_finite_difference_general(order::Int, num_stencils::Int, function_name::String)
    # Ensure num_stencils is even and valid
    if num_stencils % 2 != 0 || num_stencils < 2
        error("num_stencils must be an even number greater than or equal to 2")
    end

    # Define symmetric stencil points
    half_points = num_stencils ÷ 2
    stencil_points = collect(-half_points:half_points)
    stencil_points = stencil_points[stencil_points .!= 0]  # Exclude zero

    # Number of stencil points
    n = length(stencil_points)

    # Ensure derivative order is valid
    if order >= n
        error("The derivative order must be less than the number of stencil points.")
    end

    # Build the matrix A dynamically
    A = [stencil_points[j]^(i - 1) / factorial(i - 1) for i in 1:n, j in 1:n]

    # Target vector for the derivative
    b = zeros(n)
    b[order + 1] = 1.0  # Encode the derivative order

    # Solve for coefficients
    c = A \ b

    # Symmetry handling for odd-order derivatives (like first derivative)
    if order % 2 == 1
        # Odd-order: ensure coefficients respect asymmetry
        for i in 1:n
            if stencil_points[i] < 0
                c[i] = -c[findfirst(x -> x == -stencil_points[i], stencil_points)]
            end
        end
    else
        # Even-order: ensure coefficients respect symmetry
        c = (c .+ reverse(c)) / 2
    end

    # Generate the finite difference formula as Julia code
    terms = [string(c[i], " * g(", stencil_points[i], "*h)") for i in 1:n]
    numerator = join(terms, " + ")
    denominator = order == 0 ? "1" : "(h^$order)"
    function_code = """
    function $function_name(g, h)
        return ($numerator) / $denominator
    end
    """

    # Evaluate the function in the current module
    eval(Meta.parse(function_code))

    # Return useful results
    return stencil_points, c, function_code
end

# Generate the functions up to nth order with n strencils (defined within the input parameters as arguments 4 and 5)
function generate_all_finite_differences(max_order::Int, stencil::Int, prefix="g")
    function_names = String[]  # To store the names of generated functions

    for order in 0:max_order
        function_name = "$(prefix)$(order)_approx"
        push!(function_names, function_name)

        # Dynamically generate the finite difference function
        generate_finite_difference_general(order, stencil, function_name)
    end

    return function_names    # now contains names like "g0_approx", "g1_approx", "g2_approx", etc.
end

# Test
#function_names = generate_all_finite_differences(n_order, n_stencil)
#println(function_names)

function generate_gk_values(function_names::Vector{String}, calc_dyadics_over_beta::Function, h::Float64)
    gk_values = []  # To store the function results and preserve order

    # Loop over the function names
    for fn_name in function_names
        
        # Convert the string to a symbol and use eval to get the actual function
        fn_symbol = Symbol(fn_name)  # Convert string to symbol (function name)
        
        # Use `eval` to get the function object, then call it with the arguments
        approximation_result = eval(fn_symbol)(calc_dyadics_over_beta, h)

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

    open(joinpath(output_dir, "diff_total.txt"), "a") do f7
        println(f7, diff_norm_beta)
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

function generate_output_file()
    # Replace "." with "_" in ARGS
    args = ARGS
    arg1, arg2, arg3, arg4 = (replace(arg, "." => "_") for arg in args)

    # Generate a unique file name
    output_file = joinpath(output_dir, "resultsX$(arg1)X$(arg2)X$(arg3)X$(arg4).txt")

    # Open the file to write the output
    open(output_file, "w") do file
        # Write input parameters
        #println(file, "Input parameters:")
        #println(file, @sprintf("  S = %.3f", S))
        #println(file, @sprintf("  T = %.3f", T))
        #println(file, @sprintf("  h = %.6f", h))
        #println(file, @sprintf("  nth order = %.1f", n_order))
        #println(file, @sprintf("  number of stencils = %.1f", n_stencil))
        #println(file, "  Output File = $output_file")
        #println(file)

        # Write matrices
        println(file, "\nLucas' Dyadics:")
        write_matrix(file, "SS_beta:", SS_beta)
        write_matrix(file, "SS_beta2:", SS_beta2)
        write_matrix(file, "SS_beta3:", SS_beta3)

        println(file, "\nFull Term Derivative Dyadics:")
        write_matrix(file, "Full Term Dyadic 1:", full_beta_derivatives[1])
        write_matrix(file, "Full Term Dyadic 2:", full_beta_derivatives[2])
        write_matrix(file, "Full Term Dyadic 3:", full_beta_derivatives[3])
        write_matrix(file, "Full Term Dyadic 4:", full_beta_derivatives[4])
        write_matrix(file, "Full Term Dyadic 5:", full_beta_derivatives[5])
        write_matrix(file, "Full Term Dyadic 6:", full_beta_derivatives[6])

        println(file, "\nDyadics Comparison:")
        write_matrix(file, "Exact Dyadics:", exact_dyadics)
        write_matrix(file, "SS_Lucas:", SS_Lucas)
        write_matrix(file, "Numerical Dyadics with Full Terms:", full_beta_derivative_total)
        write_matrix(file, "Residual term (exact - SS_Lucas):", residual)
        write_matrix(file, "Numerical Residual(4th + 5th + 6th):", full_beta_residual)

        println("Results have been saved!")
    end
end

