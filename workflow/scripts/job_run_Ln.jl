using JSON
script_dir = @__DIR__
input_dir = joinpath(script_dir, "..", "data")
output_dir = joinpath(script_dir, "..", "results")

# Define a function that reads the optimal h value and uses it as an argument
function extract_h_value_from_file(filename)
    optimal_h_value = nothing  # Initialise the variable

    # Open the file and read the lines
    open(joinpath(output_dir, filename), "r") do file
        for line in eachline(file)
            # Look for the specific line containing the h value
            if occursin(r"Optimal h value \(average of all beta terms\):", line)
                # Extract the number using a regular expression
                m = match(r"[-+]?\d*\.\d+e[-+]?\d+", line)
                if m !== nothing
                    # Parse the extracted value as Float64
                    optimal_h_value = parse(Float64, m.match)
                    break  # Exit the loop 
                end
            end
        end
    end
    
    return optimal_h_value
end

function read_range(filename)
    # Read the file as a single string and remove brackets
    data = strip(read(filename, String), ['[', ']', '\n', ' '])
    
    # Split the string by commas and convert to Float64
    return parse.(Float64, strip.(split(data, ",")))
end

# Define range
Tmin = 100
Tmax = 500
Tinterval = 50

# Extract the optimal h value and assign it to arg3_range
h_value = extract_h_value_from_file(joinpath(output_dir, "optimal_h_values.txt"))
if isnothing(h_value)
    println("Failed to extract optimal h value.")
    exit()
end

# Ensuring that derivatives calculations are only done once
ENV["DERIVATIVES_DONE"] = "false"

# Build the command
command = `julia ./scripts/dyadics.jl $Tmin $Tmax $Tinterval $h_value`

println("Running: ", command)
run(command)