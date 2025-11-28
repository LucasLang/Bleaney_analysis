using MagFieldLFT
using LinearAlgebra

Ln = "Yb"
shparam = MagFieldLFT.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
sh = MagFieldLFT.SpinHamiltonian(shparam)
T = 298.0
susc = MagFieldLFT.calc_susceptibility_vanVleck(sh, T)
susc = susc-tr(susc)/3*Matrix(1.0I, 3, 3)   # traceless part of the susceptibility tensor
vals = eigvals(susc)
vecs = eigvecs(susc)
println("Ln: ", Ln)
display(vals)
display(vecs)

au2angstrom3 = (MagFieldLFT.a0/MagFieldLFT.angstrom)^3
vals *= au2angstrom3

order = sortperm(vals, by=abs)
chi_x = vals[order[1]]
chi_y = vals[order[2]]
chi_z = vals[order[3]]

chi_ax = chi_z - 0.5*(chi_x+chi_y)
chi_rh = 0.5*(chi_x - chi_y)

ez = vecs[:,order[3]]
beta = acos(ez[3]/norm(ez))
alpha = mod(atan(ez[2], ez[1]), 2pi)


rhombicity = chi_rh/chi_ax

#println("Ln: ", Ln)
#println("chi_ax (angstrom^3): ", chi_ax)
#println("rhombicity chi_rh/chi_ax: ", rhombicity)
#println("alpha: ", alpha*360/2/pi)
#println("beta: ", beta*360/2/pi)
