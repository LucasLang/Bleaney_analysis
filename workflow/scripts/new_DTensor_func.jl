using Printf
using LinearAlgebra
using DelimitedFiles

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
    temp_ratio = [3e-2, 1e-1, 3e-1, 1.0, 3.0, 10.0]  
else
    temp_ratio = [-3e-2, -1e-1, -3e-1, -1.0, -3.0, -10.0] # if D < 0
end

T = generate_temperature.(temp_ratio, D)

script_dir = @__DIR__
input_dir = joinpath(script_dir, "..", "data")
output_dir = joinpath(script_dir, "..", "results")

# Generate a unique file name
output_file = ARGS[4]
output_file2 = ARGS[5]

writedlm(output_file, Dtensor)
writedlm(output_file2, T)