using MagFieldLFT
using LinearAlgebra

Ln = "Dy"
shparam = MagFieldLFT.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
sh = MagFieldLFT.SpinHamiltonian(shparam)

T = 298.0

exact_dyadic = MagFieldLFT.calc_dyadic(sh, T)

deriv1 = MagFieldLFT.JJbeta(shparam)
deriv2 = MagFieldLFT.JJbeta2(shparam)
deriv3 = MagFieldLFT.JJbeta3(shparam)

beta = 1/MagFieldLFT.kB/T

approx1 = deriv1*beta
approx2 = approx1 + deriv2*beta^2/2
approx3 = approx2 + deriv3*beta^3/6

println("Relative error order beta: ", norm(approx1-exact_dyadic)/norm(exact_dyadic))
println("Relative error order beta^2 (Bleaney): ", norm(approx2-exact_dyadic)/norm(exact_dyadic))
println("Relative error order beta^3: ", norm(approx3-exact_dyadic)/norm(exact_dyadic))

eigvals_exact = eigvals(exact_dyadic - tr(exact_dyadic)/3*Matrix(1.0I, 3, 3))
eigvals_approx2 = eigvals(approx2 - tr(approx2)/3*Matrix(1.0I, 3, 3))
eigvals_approx3 = eigvals(approx3 - tr(approx3)/3*Matrix(1.0I, 3, 3))

println("Exact eigenvalues of anisotropic part: ", eigvals_exact)
println("Eigenvalues of anisotropic part up to order beta^2 (Bleaney's theory): ", eigvals_approx2)
println("Eigenvalues of anisotropic part up to order beta^3: ", eigvals_approx3)
