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
j1 = 1/2
m1 = -0.5
j2 = 1/2
m2 = -0.5
j = 1
m = -1
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


