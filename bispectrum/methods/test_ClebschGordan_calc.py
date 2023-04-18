"""
Compare the Clebsch Gordan Coefficient calculation from class function (calc) vs SymPy
Compare the coupling coefficient calculation from class function (calc) vs SymPy
"""
from sympy.physics.quantum.cg import CG
from sympy import *
from bispectrum.methods.calc.ClebschGordan import Clebsch_Gordan, H_coeff
import timeit
from sympy.physics.quantum.cg import CG
#Function using Sympy
j1 = 1.5
m1 = 0.5
j2 = 1
m2 = -1
j = 2
m = -0.5
t0=timeit.default_timer()
cg = CG(j1,m1,j2,m2,j,m)
cg = cg.doit()
t1=timeit.default_timer()
print(N(cg))
print("Execution time for CG function from Sympy:", t1-t0, "seconds")

#Our function
t2=timeit.default_timer()
CG_calc = Clebsch_Gordan(j1,j2,j,m1,m2,m)
cg_calc = CG_calc.cg()
print (cg_calc)
t3=timeit.default_timer()
print("Execution time for CG function from Clebsch_Gordan:", t3-t2, "seconds")
print("Execution time for CG calculation using class method is", round((t3-t2)/(t1-t0)), \
      "times faster than Sympy function")


#Coupling Coefficient Calculation
def getCoeffH(j1,j2,j,m1,m2,m,m1p,m2p,mp):
      '''
        Calculate the coupling coefficient H use SymPy
      '''
      cg = CG(j1,m1,j2,m2,j,m)
      cg = cg.doit()
      cg_p = CG(j1,m1p,j2,m2p,j,mp)
      cg_p = cg_p.doit()
      H_coeff = (cg)*(cg_p)
      H = N(H_coeff)
      return H

#Example (𝑗1,𝑗2,𝑗,𝑚1,𝑚2,𝑚,𝑚′1,𝑚′2,𝑚′)=(1,1.5,2.5,1.0,0.5,1.5,−1.0,−0.5,−1.5)
H= H_coeff(1,1.5,2.5,1.0,0.5,1.5,-1.0,-0.5,-1.5)
print ("Calc function for calculate coupling coefficient", H)

H_sympy = getCoeffH(1,1.5,2.5,1.0,0.5,1.5,-1.0,-0.5,-1.5)
print ("Sympy function for calculate coupling coefficient", H_sympy)