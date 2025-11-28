using JSON
script_dir = @__DIR__
input_dir = joinpath(script_dir, "..", "data")
output_dir = joinpath(script_dir, "..", "results")

# Function to generate runs dynamically
function generate_runs(arg1_range, arg2_range)
    runs = []
    for arg1 in arg1_range
        for arg2 in arg2_range
            
            push!(runs, (arg1, arg2, 5, 12))  # T, h, derivative order, stencil
        
        end
    end
    return runs
end

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
arg1_range = 100:50:500
#arg1_range = read_range(joinpath(input_dir, "T.txt")) 
#arg3_range = [200.0]

# Extract the optimal h value and assign it to arg3_range
h_value = extract_h_value_from_file(joinpath(output_dir, "optimal_h_values.txt"))
if isnothing(h_value)
    println("Failed to extract optimal h value.")
    exit()
end

# Use the extracted value as a single number for arg3_range
arg2_range = [h_value]

# Generate the runs
runs = generate_runs(arg1_range, arg2_range)

# Ensuring that derivatives calculations are only done once
ENV["DERIVATIVES_DONE"] = "false"

# Loop through each run
for (arg1, arg2, arg3, arg4) in runs
    # Build the command
    command = `julia ./scripts/dyadics.jl $arg1 $arg2 $arg3 $arg4`

    println("Running: ", command)
    run(command)
end

# More parameter extraction
open(joinpath(output_dir, "runs.json"), "w") do file
    write(file, JSON.json(runs))
end

# Check for txt plot points files
required_files = [
    joinpath(output_dir, "diff_norm1.txt"),
    joinpath(output_dir, "diff_norm2.txt"),
    joinpath(output_dir, "diff_norm3.txt"),
    joinpath(output_dir, "diff_norm4.txt"),
    joinpath(output_dir, "diff_norm5.txt"),
    joinpath(output_dir, "diff_norm6.txt")
]  # List of required files





