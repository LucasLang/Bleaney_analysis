using JSON

script_dir = @__DIR__
output_dir = joinpath(script_dir, "..", "results")

function generate_hruns(arg1_range, arg2_range)
    hruns = []
    for arg1 in arg1_range
        for arg2 in arg2_range
            
            push!(hruns, (arg1, arg2, 3, 12))  # T, h, derivative order, stencil
            
        end
    end
    return hruns
end

# Define range
#arg1_range = [1]
arg1_range = [0.0]
arg2_range = 10.0 .^ (-5:1:4)

# Generate the runs
hruns = generate_hruns(arg1_range, arg2_range)

# Loop through each run
for (arg1, arg2, arg3, arg4) in hruns
    # Build the command
    command = `julia ./scripts/h_sans_beta.jl $arg1 $arg2 $arg3 $arg4`

    println("Running: ", command)
    run(command)
end

# More parameter extraction
open(joinpath(output_dir, "hruns.json"), "w") do file
    write(file, JSON.json(hruns))
end

# Check for txt plot points files
required_files = [
    joinpath(output_dir, "h_norm1.txt"),
    joinpath(output_dir, "h_norm2.txt"),
    joinpath(output_dir, "h_norm3.txt")
    # joinpath(output_dir, "h_norm4.txt"),
    # joinpath(output_dir, "h_norm5.txt"),
    # joinpath(output_dir, "h_norm6.txt")
]


