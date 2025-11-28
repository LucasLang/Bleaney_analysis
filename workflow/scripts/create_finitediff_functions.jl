num_stencils = parse(Int64, ARGS[1])
outfile = ARGS[2]

max_order = num_stencils-1 
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

# Build the matrix A dynamically
A = [stencil_points[j]^(i - 1) / factorial(i - 1) for i in 1:n, j in 1:n]

function_expressions = []
function_names = String[]
file_content = ""
for order in 0:max_order
    function_name = "g$(order)_approx"
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
    global file_content *= function_code

    push!(function_names, function_name)
end
function_name_string = "function_names = [" * join(function_names, ",") * "]"
file_content *= function_name_string

write(outfile, file_content)