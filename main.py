from sympy.physics.quantum import OrthogonalKet, qapply
from sympy.physics.quantum.operatorset import operators_to_state
from sympy.physics.quantum.dagger import Dagger
from sympy import cse
from sympy import I as j
import numpy
from sympy import exp
from sympy import Symbol, simplify

from math import sqrt
zero    = OrthogonalKet(0)
one     = OrthogonalKet(1)
two     = OrthogonalKet(2)
three   = OrthogonalKet(3)
four    = OrthogonalKet(4)
five    = OrthogonalKet(5)
six     = OrthogonalKet(6)
seven   = OrthogonalKet(7)

phi = Symbol("φ", positive=True)
theta = Symbol("θ", positive=True)

def phase(angle):
    return exp(-j * angle)

def getNorm(x, y):
    return sqrt(simplify(x ** 2 + y ** 2))
    
def proj(vec):
    return vec * vec.dual

Ubs1 = 1/sqrt(2)*(j * two + three) * zero.dual + 1/sqrt(2)*(two + j * three) * one.dual
Um1 = j * four * two.dual + proj(three)
Um2 = j * five * three.dual + proj(two)

Uphi = phase(phi) * proj(three)
Utheta = phase(theta) * proj(five)

Ubs2 = 1/sqrt(2)*(j * six + seven) * four.dual + 1/sqrt(2)*(six + j * seven) * five.dual
D1 = proj(six)
D2 = proj(seven)

U = qapply(Ubs2*(Utheta*(Um2*(Uphi*(Um1*Ubs1)))))


psi = 1 / sqrt(2) * (zero + one)
print(psi)
print(simplify(qapply(U*psi)))
rho_d1 = qapply(six.dual *(D1 * (U * psi))) 
rho_d2 = qapply(seven.dual * (D2 * (U * psi) ) )
x1, y1 = rho_d1.as_real_imag()
x2, y2 = rho_d2.as_real_imag()

print(simplify(rho_d1))
print(simplify(rho_d2))

print("Probability of |6>:{}".format(simplify(getNorm(x1, y1)) * 100))
print("Probability of |7>:{}".format(simplify(getNorm(x2, y2)) * 100))


