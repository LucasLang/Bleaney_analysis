using JSON

# Global variables to store user choices
global do_norms = ""
global do_matrices = ""

# Function to prompt the user for specific choices
function prompt_user(operation_name)
    println("Would you like to run the $operation_name script? (yes/no)")
    user_choice = readline()
    while !(user_choice in ["yes", "no"])
        println("Invalid input! Please enter 'yes' or 'no'.")
        user_choice = readline()
    end
    return user_choice
end

# Initialise the user choices
function initialise_choices()
    global do_norms = prompt_user("plot points")
    global do_matrices = prompt_user("dyadics")
end

# Check if no operations are selected
function check_exit_conditions()
    if do_norms == "no" && do_matrices == "no"
        println("No operations selected. Exiting Programme...")
        exit()
    end
end

# Function to generate runs dynamically
function generate_runs(arg2_range, arg3_range)
    runs = []
    for arg1 in arg1_range
        for arg2 in arg2_range
            for arg3 in arg3_range
                push!(runs, (arg1, arg2, arg3, 5, 12))  # S, T, h, derivative order, stencil
            end
        end
    end
    return runs
end

# Define range
arg1_range = [1.0]
arg2_range = [298.0]
arg3_range = 10.0 .^(-4:1:4)

# Generate the runs
runs = generate_runs(arg2_range, arg3_range)

# Prompt user once for choices
initialise_choices()

# Check if program should exit early
check_exit_conditions()

# Loop through each run
for (arg1, arg2, arg3, arg4, arg5) in runs
    # Build the command, passing do_norms and do_matrices as environment variables
    command = `julia ./dyadics.jl $arg1 $arg2 $arg3 $arg4 $arg5`
    
    println("Running: ", command)
    println("Environment: Plots=$do_norms, Dyadics=$do_matrices")

    # Execute the command with the environment variables set
    withenv("DO_NORMS" => do_norms, "DO_MATRICES" => do_matrices) do
        run(command)
    end
end

# More parameter extraction
open("runs.json", "w") do file
    write(file, JSON.json(runs))
end

# Check for txt plot points files
required_files = ["diff_norm1.txt", "diff_norm2.txt", 
                  "diff_norm3.txt", "diff_norm4.txt", 
                  "diff_norm5.txt", "diff_norm6.txt"]  

if all(isfile, required_files)
    println("Plot points found. Generating error plot...")
    
    # Run the Python script
    run(`py h2_plots.py`)
    println("Error plot generated!")
else
    missing_files = filter(f -> !isfile(f), required_files)
    println("The following required file(s) are missing: ", missing_files)
    println("Diagram will not be generated.")
end

# List of files to delete
files_to_delete = ["diff_norm1.txt", "diff_norm2.txt", 
                  "diff_norm3.txt", "diff_norm4.txt", 
                  "diff_norm5.txt", "diff_norm6.txt",
                  "diff_total.txt", "runs.json"]

# Functions to delete files
function cleanup_files()
    for file in files_to_delete
        if isfile(file)
            rm(file)
            println("Deleted: $file")
        else
            println("File not found: $file")
        end
    end
end

function end_clean_prompt()
    println("Do you want to clean up the input files? (yes/no)")
    response = readline()
    if response == "yes"
        cleanup_files()
    else
        println("Plot iput files retained.")
    end
end

atexit(end_clean_prompt)
