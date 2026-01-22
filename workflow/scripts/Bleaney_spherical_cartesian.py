import sympy as sp
from sympy.physics.quantum.cg import CG

Dxx, Dyy, Dxy, Dxz, Dyz = sp.symbols("Dxx Dyy Dxy Dxz Dyz")
Dzz = -Dxx-Dyy

Dtensor = sp.Matrix([[Dxx, Dxy, Dxz], [Dxy, Dyy, Dyz], [Dxz, Dyz, Dzz]])   # Dtensor is symmetric and traceless

# matrix A for conversion between Cartesian and spherical components of a vector
A = sp.Matrix([[-1/sp.sqrt(2), -sp.I/sp.sqrt(2), 0], [0, 0, 1], [1/sp.sqrt(2), -sp.I/sp.sqrt(2), 0]])  
Ncomp_vector = 3

ind = {+1:0, 0:1, -1:2}  # map "magnetic quantum numbers" to indices 0,1,2 of the SymPy matrix

# this function takes a 3x3 matrix (2nd rank tensor) and returns the corresponding rank-2 spherical tensor component Q
def calc_Cart2Spher_tensor(Dtensor, Q):
    Duncoupled = A*Dtensor*A.T
    result = 0
    for q in range(-1, -1+Ncomp_vector):
        for qprime in range(-1, -1+Ncomp_vector):
            result += CG(1, q, 1, qprime, 2, Q).doit() * Duncoupled[ind[q],ind[qprime]]
    return sp.simplify(result)

Ncomp_rank2tensor = 2*2+1   # multiplicity of a second-order spherical tensor = 2L+1 = 5

D2 = {Q:calc_Cart2Spher_tensor(Dtensor, Q) for Q in range(-2, -2+Ncomp_rank2tensor)}

def calc_D2xD2(K, Q):
    result = 0
    for q in range(-2, -2+Ncomp_rank2tensor):
        for qprime in range(-2, -2+Ncomp_rank2tensor):
            result += CG(2, q, 2, qprime, K, Q).doit() * D2[q] * D2[qprime]
    return result


D2D2_2 = {Q: calc_D2xD2(2, Q) for Q in range(-2,-2+Ncomp_rank2tensor)}

D2D2_2_uncoupled = sp.zeros(3,3)
for q in range(-1, -1+Ncomp_vector):
    for qprime in range(-1, -1+Ncomp_vector):
        for Q in range(-2, -2+Ncomp_rank2tensor):
            D2D2_2_uncoupled[ind[q], ind[qprime]] += CG(1, q, 1, qprime, 2, Q).doit() * D2D2_2[Q]

D2D2_2_Cartesian = A.H*D2D2_2_uncoupled*A.T.H

Dsquared = Dtensor * Dtensor

Dsquared_aniso = Dsquared - (Dsquared.trace()/3)*sp.Matrix([[1,0,0], [0,1,0], [0,0,1]])

print("\nRatio of the Cartesian form of (D2xD2)^2 and (D^2)^aniso:")
for i in range(0,3):
    for j in range(i,3):
        print(f"Element ({i},{j}): ", sp.simplify(sp.expand(D2D2_2_Cartesian[i,j])/ sp.expand(Dsquared_aniso[i,j])))

D2D2_0 = calc_D2xD2(0, 0)
B2B2_0_uncoupled = sp.zeros(3,3)
ind = {+1:0, 0:1, -1:2}
for q in range(-1, -1+Ncomp_vector):
    for qprime in range(-1, -1+Ncomp_vector):
        B2B2_0_uncoupled[ind[q], ind[qprime]] += CG(1, q, 1, qprime, 0, 0).doit() * D2D2_0

B2B2_0_Cartesian = A.H*B2B2_0_uncoupled*A.T.H

tr_Dsquared_times_I = Dsquared.trace()*sp.Matrix([[1,0,0], [0,1,0], [0,0,1]])

print("\nRatio of the Cartesian form of (D2xD2)^0 and tr(D^2)*I (for diagonal elements):")
for i in range(0,3):
    print(f"Element ({i},{i}): ", sp.simplify(sp.expand(B2B2_0_Cartesian[i,i]) / tr_Dsquared_times_I[i,i]))

print("\nOff-diagonal elements of the Cartesian form of (D2xD2)^0:")
for i in range(0,3):
    for j in range(i+1,3):
        print(f"Element ({i},{j}): ", sp.simplify(sp.expand(B2B2_0_Cartesian[i,j])))