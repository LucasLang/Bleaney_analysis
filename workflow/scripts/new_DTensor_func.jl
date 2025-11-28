using Printf
using LinearAlgebra
using JSON

function new_Dtensor(D, E)
    conversion_factor = 1 / 219474.63  # cm^-1 to a.u. conversion
    return D * conversion_factor * [
        (-1/3 + E/D)   0             0;
        0            (-1/3 - E/D)    0;
        0             0             2/3
    ]
end

D = parse(Float64, ARGS[1])  
E = parse(Float64, ARGS[2])
S = parse(Float64, ARGS[3])

Dtensor = new_Dtensor(D, E)

ratio = E/D
k = 3.166811563e-6 
conversion_factor_k = 219474.63
kB = k * conversion_factor_k

# temp_ratio = D/kB*T

function generate_temperature(temp_ratio, D)
    conversion_factor = 1 / 219474.63
    result = (D * conversion_factor) / (k * temp_ratio)
    return abs(result)  # Ensuring that T is positive
end

if D > 0
    temp_ratio = [3e-2, 1e-1, 3e-1, 1.0, 3.0, 10.0] #10.0 .^ (-2:1:1) 
else
    temp_ratio = [-3e-2, -1e-1, -3e-1, -1.0, -3.0, -10.0] #-10.0 .^ (-2:1:1) # if D < 0
end

T = generate_temperature.(temp_ratio, D)

script_dir = @__DIR__
input_dir = joinpath(script_dir, "..", "data")
output_dir = joinpath(script_dir, "..", "results")

function write_matrix(file, matrix)
    #println(file, "\n$name:")
    for i in 1:size(matrix, 1)
        for j in 1:size(matrix, 2)
            @printf(file, " %12.6e ", matrix[i, j])
        end
        println(file)
    end
end

function generate_Dtensor_file()
    # Replace "." with "_" in ARGS
    args = ARGS
    arg1, arg2 = (replace(arg, "." => "_") for arg in args)

    # Generate a unique file name
    output_file = joinpath(input_dir, "spin_dyadics_input.txt")
    output_file2 = joinpath(input_dir, "S.txt")
    output_file3 = joinpath(input_dir, "D.txt")
    output_file4 = joinpath(input_dir, "T.txt")
    output_file5 = joinpath(input_dir, "ratio.txt")

    open(output_file, "w") do file
        write_matrix(file, Dtensor)

    end
   
    open(output_file2, "w") do file
        println(file, S)

    end

    open(output_file3, "w") do file
        println(file, D)

    end
    
    open(output_file4, "w") do file
        println(file, T)

    end
    
    open(output_file5, "w") do file
        println(file, ratio)

    end

end

generate_Dtensor_file()