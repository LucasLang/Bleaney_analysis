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

# Build the command
command = `julia ./scripts/h_sans_beta.jl`

println("Running: ", command)
run(command)

# More parameter extraction
open(joinpath(output_dir, "hruns.json"), "w") do file
    write(file, JSON.json(hruns))
end