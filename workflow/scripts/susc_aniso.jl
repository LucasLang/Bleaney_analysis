using ParaMag
using LinearAlgebra

datafolder = ARGS[1]
outfile = ARGS[2]
Ln_list = ARGS[3:end]

atomic_numbers = Dict("La" => 57, "Ce" => 58, "Pr" => 59, "Nd" => 60, "Pm" => 61,
                      "Sm" => 62, "Eu" => 63, "Gd" => 64, "Tb" => 65, "Dy" => 66,
                      "Ho" => 67, "Er" => 68, "Tm" => 69, "Yb" => 70)

Ln_list_sorted = sort(Ln_list, by = x -> atomic_numbers[x])

function calc_chi_ax_rhombicity_fromsusc(susc)
    susc_traceless = susc-tr(susc)/3*Matrix(1.0I, 3, 3)
    vals = eigvals(susc_traceless)
    vecs = eigvecs(susc_traceless)

    au2angstrom3 = (ParaMag.a0/ParaMag.angstrom)^3
    vals *= au2angstrom3

    order = sortperm(vals, by=abs)
    chi_x = vals[order[1]]
    chi_y = vals[order[2]]
    chi_z = vals[order[3]]

    chi_ax = chi_z - 0.5*(chi_x+chi_y)
    chi_rh = 0.5*(chi_x - chi_y)

    rhombicity = chi_rh/chi_ax
    return chi_ax, rhombicity
end

function chi_ax_rhombicity_Bleaney(shparam, T)
    dyadic = ParaMag.calc_dyadic_order2(shparam, T)
    susc = ParaMag.calc_susceptibility_fromdyadic(dyadic, shparam.gtensor)
    return calc_chi_ax_rhombicity_fromsusc(susc)
end

function chi_ax_rhombicity_3rdorder(shparam, T)
    dyadic = ParaMag.calc_dyadic_order3(shparam, T)
    susc = ParaMag.calc_susceptibility_fromdyadic(dyadic, shparam.gtensor)
    return calc_chi_ax_rhombicity_fromsusc(susc)
end

function chi_ax_rhombicity_exact(sh, T)
    susc = ParaMag.calc_susceptibility_vanVleck(sh, T)
    return calc_chi_ax_rhombicity_fromsusc(susc)
end

Ln_number = length(Ln_list_sorted)

# Column 1: Bleaney. Column 2: 3rd order. Column 3: exact (infinite order).
chi_ax_matrix = Array{Float64}(undef, Ln_number, 3)
chi_rh_matrix = Array{Float64}(undef, Ln_number, 3)

T = 298.0

for (i,Ln) in enumerate(Ln_list_sorted)
    shparam = ParaMag.SHParam_lanthanoid("$(datafolder)/Bkq_$(Ln)_real", Ln)
    sh = ParaMag.SpinHamiltonian(shparam)
    chi_ax_matrix[i,1], chi_rh_matrix[i,1] = chi_ax_rhombicity_Bleaney(shparam, T)
    chi_ax_matrix[i,2], chi_rh_matrix[i,2] = chi_ax_rhombicity_3rdorder(shparam, T)
    chi_ax_matrix[i,3], chi_rh_matrix[i,3] = chi_ax_rhombicity_exact(sh, T)
end

println(Ln_list_sorted)
display(chi_ax_matrix)
display(chi_rh_matrix)
