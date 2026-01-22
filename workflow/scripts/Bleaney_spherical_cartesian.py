import sympy as sp
from sympy.physics.quantum.cg import CG

Dxx, Dyy, Dxy, Dxz, Dyz = sp.symbols("Dxx Dyy Dxy Dxz Dyz")
Dzz = -Dxx-Dyy

Dtensor = sp.Matrix([[Dxx, Dxy, Dxz], [Dxy, Dyy, Dyz], [Dxz, Dyz, Dzz]])   # Dtensor is symmetric and traceless

A = sp.Matrix([[-1/sp.sqrt(2), -sp.I/sp.sqrt(2), 0], [0, 0, 1], [1/sp.sqrt(2), -sp.I/sp.sqrt(2), 0]])
mult_triplet = 3

def calc_Cart2Spher_tensor(Dtensor, Q):
    Duncoupled = A*Dtensor*A.T
    result = 0
    ind = {+1:0, 0:1, -1:2}
    for q in range(-1, -1+mult_triplet):
        for qprime in range(-1, -1+mult_triplet):
            result += CG(1, q, 1, qprime, 2, Q).doit() * Duncoupled[ind[q],ind[qprime]]
    return sp.simplify(result)

#import pdb; pdb.set_trace()
#calc_Cart2Spher_tensor(Dtensor, +2)




# B2_2 = (Dxx-Dyy+ 2*sp.I*Dxy)/2
# B2_1 = -Dxz-sp.I*Dyz
# B2_0 = (2*Dzz - Dxx - Dyy)/sp.sqrt(6)
# B2_m1 = Dxz -sp.I*Dyz
# B2_m2 = (Dxx - Dyy -2*sp.I*Dxy)/2

mult_quintet = 2*2+1   # multiplicity of a second-order spherical tensor = 2L+1 = 5

B2 = {Q:calc_Cart2Spher_tensor(Dtensor, Q) for Q in range(-2, -2+mult_quintet)}
print(B2)

def calc_B2B2(K, Q):
    result = 0
    for q in range(-2, -2+mult_quintet):
        for qprime in range(-2, -2+mult_quintet):
            result += CG(2, q, 2, qprime, K, Q).doit() * B2[q] * B2[qprime]
    return result


B2B2_2 = {Q: calc_B2B2(2, Q) for Q in range(-2,-2+mult_quintet)}

# B2B2_2_xx = sp.factor(sp.simplify((B2B2_2[2] + B2B2_2[-2])/2 - B2B2_2[0]/sp.sqrt(6)))
# B2B2_2_yy = sp.factor(sp.simplify(-(B2B2_2[2] + B2B2_2[-2])/2 - B2B2_2[0]/sp.sqrt(6)))
# B2B2_2_zz = sp.factor(sp.simplify(2*B2B2_2[0]/sp.sqrt(6)))
# B2B2_2_xy = sp.factor(sp.simplify(-sp.I*(B2B2_2[2] - B2B2_2[-2])/2))
# B2B2_2_xz = sp.factor(sp.simplify(-(B2B2_2[1] - B2B2_2[-1])/2))
# B2B2_2_yz = sp.factor(sp.simplify(sp.I*(B2B2_2[1] + B2B2_2[-1])/2))

B2B2_2_uncoupled = sp.zeros(3,3)
ind = {+1:0, 0:1, -1:2}
for q in range(-1, -1+mult_triplet):
    for qprime in range(-1, -1+mult_triplet):
        for Q in range(-2, -2+mult_quintet):
            B2B2_2_uncoupled[ind[q], ind[qprime]] += CG(1, q, 1, qprime, 2, Q).doit() * B2B2_2[Q]

B2B2_2_Cartesian = A.H*B2B2_2_uncoupled*A.T.H

# Test where I just back-transform B2:
#B2B2_2_xx = sp.simplify((B2[2] + B2[-2])/2 - B2[0]/sp.sqrt(6))
#B2B2_2_yy = sp.simplify(-(B2[2] + B2[-2])/2 - B2[0]/sp.sqrt(6))
#B2B2_2_zz = sp.simplify(2*B2[0]/sp.sqrt(6))
#B2B2_2_xy = sp.simplify(-sp.I*(B2[2] - B2[-2])/2)
#B2B2_2_xz = sp.simplify(-(B2[1] - B2[-1])/2)
#B2B2_2_yz = sp.simplify(sp.I*(B2[1] + B2[-1])/2)

print("\nCartesian form of (D2xD2)^2 tensor:")
print("xx: ", sp.simplify(B2B2_2_Cartesian[0,0]))
print("yy: ", sp.simplify(B2B2_2_Cartesian[1,1]))
print("zz: ", sp.simplify(B2B2_2_Cartesian[2,2]))
print("xy: ", sp.simplify(B2B2_2_Cartesian[0,1]))
print("xz: ", sp.simplify(B2B2_2_Cartesian[0,2]))
print("yz: ", sp.simplify(B2B2_2_Cartesian[1,2]))

D2 = Dtensor * Dtensor

D2_aniso = D2 - (D2.trace()/3)*sp.Matrix([[1,0,0], [0,1,0], [0,0,1]])
print("D2_aniso xx: ", sp.simplify(sp.expand(D2_aniso[0,0])))
print("D2_aniso yy: ", sp.simplify(sp.expand(D2_aniso[1,1])))
print("D2_aniso zz: ", sp.simplify(sp.expand(D2_aniso[2,2])))
print("D2_aniso xy: ", sp.simplify(sp.expand(D2_aniso[0,1])))
print("D2_aniso xz: ", sp.simplify(sp.expand(D2_aniso[0,2])))
print("D2_aniso yz: ", sp.simplify(sp.expand(D2_aniso[1,2])))

print("Ratio xx: ", sp.simplify(sp.simplify(sp.expand(D2_aniso[0,0])) / sp.simplify(B2B2_2_Cartesian[0,0])))
print("Ratio yy: ", sp.simplify(sp.simplify(sp.expand(D2_aniso[1,1])) / sp.simplify(B2B2_2_Cartesian[1,1])))
print("Ratio zz: ", sp.simplify(sp.simplify(sp.expand(D2_aniso[2,2])) / sp.simplify(B2B2_2_Cartesian[2,2])))
print("Ratio xy: ", sp.simplify(sp.simplify(sp.expand(D2_aniso[0,1])) / sp.simplify(B2B2_2_Cartesian[0,1])))
print("Ratio xz: ", sp.simplify(sp.simplify(sp.expand(D2_aniso[0,2])) / sp.simplify(B2B2_2_Cartesian[0,2])))
print("Ratio yz: ", sp.simplify(sp.simplify(sp.expand(D2_aniso[1,2])) / sp.simplify(B2B2_2_Cartesian[1,2])))


B2B2_0 = calc_B2B2(0, 0)
B2B2_0_uncoupled = sp.zeros(3,3)
ind = {+1:0, 0:1, -1:2}
for q in range(-1, -1+mult_triplet):
    for qprime in range(-1, -1+mult_triplet):
        B2B2_0_uncoupled[ind[q], ind[qprime]] += CG(1, q, 1, qprime, 0, 0).doit() * B2B2_0

B2B2_0_Cartesian = A.H*B2B2_0_uncoupled*A.T.H
D2_iso = (D2.trace()/3)*sp.Matrix([[1,0,0], [0,1,0], [0,0,1]])

print("\nCartesian form of (D2xD2)^0 tensor:")
print("xx: ", sp.simplify(B2B2_0_Cartesian[0,0]))
print("yy: ", sp.simplify(B2B2_0_Cartesian[1,1]))
print("zz: ", sp.simplify(B2B2_0_Cartesian[2,2]))
print("xy: ", sp.simplify(B2B2_0_Cartesian[0,1]))
print("xz: ", sp.simplify(B2B2_0_Cartesian[0,2]))
print("yz: ", sp.simplify(B2B2_0_Cartesian[1,2]))

print("trace(D2): ", sp.expand(D2.trace()))
# print("D2_iso xx: ", sp.simplify(sp.expand(D2_iso[0,0])))
# print("D2_iso yy: ", sp.simplify(sp.expand(D2_iso[1,1])))
# print("D2_iso zz: ", sp.simplify(sp.expand(D2_iso[2,2])))
# print("D2_iso xy: ", sp.simplify(sp.expand(D2_iso[0,1])))
# print("D2_iso xz: ", sp.simplify(sp.expand(D2_iso[0,2])))
# print("D2_iso yz: ", sp.simplify(sp.expand(D2_iso[1,2])))
