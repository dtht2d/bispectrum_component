# Compute Bispectrum Components

**Color code**

**Definition**

**Code/Tools**

**Results/In-progress/Example**

********Idea********

**Limitations**

**Note**

---

# Table of content

### Motivation

### Project outline

### **I. Background**

**1. Local atomic density neighbors function in 3D/ 4D dimensional space** 

1.1. 3D atomic neighbor density function

1.2. 4D atomic neighbor density function around a central atom $i$ at location $**r**$

### **I**I. Methods

**1. Assumptions**

**2. Computational steps**

**2.1. Data prep and parameter**

**a. Data prep**

Read CIF file as data frame

**b. Define position of neighbor atoms $k$ relative to a central atom $i$  is a point within the 3D ball of radius**  

Estimate list of potentially atoms in the center cell

Choose a center atom $i$  from the list and get its coordinate 
$(x_i,y_i,z_i)$

Re-calculate coordinate for all atoms in the cell w.r.t new reference of frame  (origin at $(x_i,y_i,z_i)$)

Compute neighbors list in a chosen cutoff radius $R_{cut}$ with respect to new reference frame where its origin is at center atom $i$  location 

**c. Map possible neighbor atoms on to the set of points $(\theta_0, \theta,\phi)$** 

Compute $\theta_0, \ \theta, \ \phi$

**2.2. Compute bispectrum component for input values $(j_1, j_2, j)$**

1. Compute Clebsch- Gordan coefficient $C^{jm}_{{j_1}{m_1}{j_2}{m_2}}, C^{jm'}_{{j_1}{m'_1}{j_2}{m'_2}}$ 
2. Compute coupling coefficient   $H^{jmm'}_{{{j_1}{m_1}{m'_1}} ,{{j_2}{m_2}{m'_2}}}$ 
3. **Compute $u^j_{mm'}, \mu_{m_1m_1'}^{j_1}, \mu_{m_2,m_2'}^{j_2}$ density coefficient**
    - Compute $U_{m,m'}^{j}(\theta_{0},\theta, \phi) \equiv U^j_{m,m'}(\omega,\Theta, \Phi)$ through Wigner -D function
    - Compute $f_c(r_{ik},R_{cut}^{ik})$ cut off function

### References

# Motivation

---

- Write Python scripts to automatically generate bispectrum coefficients from a given crystal structure dataset

Why bispectrum component? 

Ref. [3] bispectrum is a three-point correlation function, is a much richer system of invariants and can provide an almost one to one representation of the atomic neighborhood

→ Someway to quantify local environment (uniform, systematic, that works for all kind of system).

# Project outlines

---

# I. Background

**Mathematical equation for bispectrum component**

$$
B_(j_1,j_2,j) = \displaystyle\sum_{{m_1,m'_1 =-j_1}}^{j_1} \displaystyle\sum_{{m_2,m'_2 = -j_2}}^{j_2} \displaystyle\sum_{{m,m'= -j}}^j  \left(u^{\smash{j}}_{m,m'}\right)^*H^{jmm'}_{{{j_1}{m_1}{m'_1}} ,{{j_2}{m_2}{m'_2}}}  u^{j_1}_{{m_1},{m_2}} u^{j_2}_{{m_2},{m_2}}
$$

## 1. **Local atomic density neighbors function in 3D/ 4D dimensional space**

### 1**.1.** **3D atomic neighbor density function**

---

The total density of neighbor atoms around a central atom $i$  of element $\mu_i$ located at the origin is represented as:

 Ref.[1] equation, notion of Ref. [2] p.6

where 

- $\delta(r_{ik})$ is the Dirac delta function centered at each neighboring site $k$ [2]
- $f_c$  is the cutoff function to ensure smooth decay for the neighbor atomic density to zero cutoff radius $R_{cut}$ [2]
    
    $$
    \begin{align*} \rho_i\left(r\right) & = \ \delta\left(r\right) + \displaystyle\sum_{r_{ik}<R_{cut}} f_c\left(r_{ik}\right)\omega_k\delta\left(r-r_{ik}\right) \\  \rho_i\left(r\right) & = \ w_{\mu_i \mu}^{self} \ \delta(0) + \displaystyle\sum_{r_{ik}<R_{cut}^{\mu_i \mu_k}} f_c\left(r_{ik};R^{\mu_i \mu_k}_{cut}\right)\omega_{\mu_k}\delta\left(r_{ik}\right)\end{align*} \tag{0}
    $$
    
- $\omega_{\mu_k}$ are the coefficients that are dimensionless weights that are chosen to distinguish atoms of different types, while the central atom is arbitrarily assigned a unit weight. [1]
- $w_{\mu_i\mu}^{self}$ is self contribution ( either 1 or 0)

### 1.2. **4D atomic neighbor density function around a central atom $i$ at location $r$**

---

$$
\rho_i(r) = \displaystyle\sum_{j=0}^{\infty} \displaystyle\sum_{m,m'=-j}^{j}u_{m,m'}^j.U^j_{m,m'}(\theta_0,\theta,\phi)\tag{1}
$$

where $u^j_{mm'}$  and $U^j_{mm'}$ are rank $(2j+1)$ complex square matrices. The . symbol indicates the scalar product of the two matrices

The expansion density coefficients (Fourier coefficient) are given by inner product of neighbor density with the basis functions $U^j_{mm'}$ because the neighbor density is a weighted sum of $\delta$- functions, each expansion coefficient can be written as sum over discrete values of the corresponding basis function [8]

$$
u_{m,m'}^j = w^{self}_{\mu_{i}}U_{m,m'}^j(0,0,0) + \displaystyle\sum_{r_{ik}<R^{\mu_i\mu_k}_{cut}} \delta_{\mu \mu_k} f_c(r_{ik};R^{\mu_i\mu_k}_{cut})w_{\mu_{k}}U_{m,m'}^{j}(\theta_{0},\theta, \phi) \tag{2}
$$

where $\delta_{\mu \mu_k}$indicates only neighbor atom of element $\mu$ contribute to partial density $\rho^\mu$( explain more about a explicit representation of different chemical elements in method session)  

---

# II. Methods

# 1. **Assumptions**

---

Assumes a linear relationship between atom energy and bispectrum components.” (Thompson et al., p. 1)

The calculation of the bispectrum components and the SNAP potential are implemented in the LAMMPS parallel molecular dynamics code. (Thompson et al., p. 2)

Bispectrum equation define as: 

$$
B_(j_1,j_2,j) = \displaystyle\sum_{{m_1,m'_1 =-j_1}}^{j_1} \displaystyle\sum_{{m_2,m'_2 = -j_2}}^{j_2} \displaystyle\sum_{{m,m'= -j}}^j  \left(u^{\smash{j}}_{m,m'}\right)^*C^{jm}_{{j_1}{m_1}{j_2}{m_2}} \times C^{jm'}_{{j_1}{m'_1}{j_2}{m'_2}} u^{j_1}_{{m_1},{m_1'}} u^{j_2}_{{m_2},{m_2'}} \tag{1}
$$

$C^{jm}_{{j_1}{m_1}{j_2}{m_2}} C^{jm'}_{{j_1}{m'_1}{j_2}{m'_2}} = H^{jmm'}_{{{j_1}{m_1}{m'_1}} ,{{j_2}{m_2}{m'_2}}}$   are coupling coefficients, analogous to the Clebsch- Gordan coefficients for rotations on 2-sphere.

# 2. Computational steps

---

### Import libraries

```python
#file directory and format
import os
from pymatgen.io.cif import CifParser 
import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_d
from sympy.physics.quantum.spin import Rotation
from sympy import *
import itertools
```

## 2.1. Data prep and parameter

---

## a. Data prep

**→ set of data, neighbor list, it could be x, y, z coordinate and convert to r, $\theta_0,\theta,\phi$**

- read CIF file as data frame
    
    [avgBL-Model.cif](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/avgBL-Model.cif)
    
- cell length (unit angstrom)
- $R_{cut }$ cut off radius
    
    ********Note:******** Chosen $R_{cut}$ needs to divide by the true cell length (i.e: with example dataset we need to divide by 18.7337) since $(x,y,z)$  are fraction with cell dimension (1,1,1), $R_i, R_k$ (center, neighbor): atom radius, measures from center of nucleus to the outermost isolated electron also so need to divide by 
    
- $r_\mu$  atomic radius for all chemical elements
- $x, y, z$ coordinate for all atoms
- atom type

![Figure 0. C**rystalline unit cell of silicon, Avogradro visualization**](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/Screen_Shot_2022-10-19_at_10.12.34_PM.png)

Figure 0. C**rystalline unit cell of silicon, Avogradro visualization**

```python
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
path = "/Users/DuongHoang/UMKC-Grad/UMKC_Research/bispectrump_component/data/avgBL-Model.cif"
dico = MMCIF2Dict(path)
df_cif = pd.DataFrame.from_dict(dico, orient='index')
x = df_cif.iloc[-3]
y = df_cif.iloc[-2]
z = df_cif.iloc[-1]
atom_type = df_cif.iloc[-4]
x_array = np.array(x[0],dtype=float)
y_array = np.array(y[0],dtype=float)
z_array = np.array(z[0],dtype=float)
atom_type_array = np.array(atom_type[0], dtype=str)
df = pd.DataFrame({"atom_type":atom_type_array,"X" : x_array, "Y":y_array, "Z": z_array})
print(df)
```

**Comment**

- Define index as orientation, column label in CIF file now become row,  index 0 column store all value atom coordinate values that need to unpack and assigned as new numpy array to create new data frame that column name is X, Y, Z coordinates and row are coordinate values for each atom in the cell
    
    **Note:** $x, y, z$  are fraction of cell length to 1. True value needs to multiply with actual cell length
    

  Out [ ]

![Output 1. Data frame of input CIF file](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/Screen_Shot_2022-10-18_at_2.51.50_PM.png)

Output 1. Data frame of input CIF file

![Output 2. Date frame with atoms (x,y,z)  coordinates and atom_type](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/Screen_Shot_2022-10-27_at_6.05.59_PM.png)

Output 2. Date frame with atoms (x,y,z)  coordinates and atom_type

---

## b. **Define position of a neighbor atom $k$ relative to a central atom $i$  is a point within the 3D ball of radius**

**Idea :** choose a central atom  $i$ where its coordinate $(x_i,y_i,z_i)$ in the center of a single cell (reference coordinate frame, origin point $O(0,0,0)$ and orthogonal axes $X,Y,Z$). Re-calculate coordinate for all atoms in the cell by shift origin to $(x_i,y_i,z_i)$. 

Coordinate of neighbor atom with respect to the new system is:

$$
\begin{bmatrix}x_k\\y_k\\z_k \end{bmatrix} =\begin{bmatrix}x -x_i\\y-y_i\\z-z_i \end{bmatrix} \tag{2}
$$

![Figure 1. Coordinate transformation, shift origin O to center atom i location ](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/IMG_2046.jpg)

Figure 1. Coordinate transformation, shift origin O to center atom i location 

![Figure 2. Draw a sphere around center atom with radius cutoff](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/IMG_2047.jpg)

Figure 2. Draw a sphere around center atom with radius cutoff

**Estimate list of potentially atoms in the center cell**

**Condition:** $\text{range}_{min} + R_{cut}, \text{range}_{max}+R_{cut} \leq$  1 

****************Example:****************

$$
0.5\leq x,y,z \leq 0.7
$$

```python
#Estimate list of potentially atoms in the center cell
df_atoms = df[(df['X'].between(0.5,0.7,inclusive='both'))
                         & (df['Y'].between(0.5,0.7,inclusive='both'))
                         & (df['Z'].between(0.5,0.7, inclusive='both'))]
print (df_atoms)
```

**Note:** the location ranges for $x, y, z$ can change in the case we only want to deal with atoms, in the range that when applying a spherical space with $R_{cut}$ around them, the sphere would not goes over the  edge’s cell. For those atoms  near the edge’s cell, different location calculation  need to apply.

![Output  3. Atoms in center cell in range x,y,z  =[0.5,0.7]](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/Screen_Shot_2022-10-27_at_6.40.03_PM.png)

Output  3. Atoms in center cell in range x,y,z  =[0.5,0.7]

**Choose a center atom $i$  from the list and get its coordinate 
$(x_i,y_i,z_i)$**

```python
#Choose a center atom i, in this example we choose atom 'Name'=17 from df_atoms dataframe
atom_i =df.iloc[17]
print (atom_i)
# id
x_i = df['X'].iloc[17]
y_i = df['Y'].iloc[17]
z_i = df['Z'].iloc[17]
print(x_i,y_i,z_i)
```

**Note:** Choose atom ID 17 as center atom, in df_atoms its index  row is [0] but in df_new its index row is [17]. 

**Re-calculate coordinate for all atoms (w.r.t new reference of frame  (origin at $(x_i,y_i,z_i)$), and distance $r_{ik}$  from a center atom $i$  to other atom in the cell**

```python
X_array = df['X'].to_numpy()
Y_array = df['Y'].to_numpy()
Z_array = df['Z'].to_numpy()
X_k_array = X_array - x_i
Y_k_array = Y_array - y_i
Z_k_array = Z_array - z_i
r_ik= np.sqrt(np.square(X_k_array)+np.square(Y_k_array)+np.square(Z_k_array))
df['X_k'],df['Y_k'], df['Z_k'],df['r_ik']= X_k_array, Y_k_array, Z_k_array, r_ik
print(df)
#Check to see if chosen center atom coordinate sets to (0,0,0)
print(df.iloc[17])
```

Check to see if chosen center atom coordinate sets to (0,0,0) 

**Compute neighbors list in a chosen cutoff radius $R_{cut}$ with respect to new reference frame where its origin is at center atom $i$  location** 

**Our approach:** 

$$
R_{cut}= some \ reasonable \ value
$$

**Example:** Sillicon atomic radius is $1.46\text{\AA}$ , Silicon–silicon π single bond ****$2.853 \text{\AA}$**, $R_{cut}= \frac{1.46+1.46}{cell\ length}=\frac{2.92}{18.7337} \approx 0.16$ but let choose $R_{cut} =0.25$

$$
r_{ik}= \sqrt{x_k^2 +y_k^2+z_k^2} \tag{3}
$$

Assume atom of element $\mu$ (in this case it’s Silicon) is a sphere with radius $r_\mu$. We want to make sure that for the neighbors list only atoms whose entire sphere is in within the bigger sphere with radius $R_{cut}$. 

$$
r_\mu = \frac{atomic\ radius \ (\text{\AA})}{cell \ length \ (\text{\AA})} = \frac{1.46}{18.7337} \approx 0.0779
$$

**Condition:**

$$
r_{ik}+r_\mu\leq R_{cut} \tag{4}
$$

```python
#INPUT values
atomic_radius = 1.46            #silicon atomic radius, unit: angstrom
cell_length = df.iloc[4]        #index row start from 0 _cell_length_a at row 5 index [4]

r_mu = 0.0779                   #scale atomic radius w.r.t cell length
R_cut = 0.25                     #scaled value w.r.t cell length (for Si-Si case)
df_ik = df[(df['r_ik'] + r_mu)<= (R_cut)].copy(deep=true)
print(df_ik[['X_k', 'Y_k', 'Z_k', 'r_ik']])
```

LAMMPS approach

$$
R_{ik}= r_{cutfrac}(R_i+R_k)
$$

**Note:** LAMMPS add $r_{cutfrac} =$  scale factor applied to all cutoff radii (positive real) (i.e. )

$$
R_{cut}=R_{ik} =(R_i+R_k)
$$

![Figure 3. Visualize cut off region around center atom i](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/IMG_1991.jpg)

Figure 3. Visualize cut off region around center atom i

![Output 4.1. Neighbors atom coordinate in the center atom’s frame of reference, and distance from the center atom to its neighbors in the cut off region.](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/Screen_Shot_2022-10-26_at_3.51.39_PM.png)

Output 4.1. Neighbors atom coordinate in the center atom’s frame of reference, and distance from the center atom to its neighbors in the cut off region.

Let try to locate  center atom to see if the  neighbors list makes sense

![Output 4.2 Center and neighbors atom (x,y,z) coordinate in cell’s frame of reference   ](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/Screen_Shot_2022-11-02_at_8.46.54_AM.png)

Output 4.2 Center and neighbors atom (x,y,z) coordinate in cell’s frame of reference   

![IMG_2055.jpg](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/IMG_2055.jpg)

Question: 

- why we missing atom id310?
- Do we need to adjust $R_{cut}$? What would be a reasonable value?
- What does it mean when we included neighbor atom that center atom does have bond with? Can we predict crystal structure  from bispectrum component?

---

## c. Map possible neighbor atoms on to the set of points $(\theta_0, \theta,\phi)$

Ref.[3], eq.(3), p. 2

**Idea:** project the atomic density onto the surface of four- dimensional unit sphere, similarly to how the Riemann sphere is constructed with the transformation

$$
(\theta_0, \theta,\phi) =[|r|/r_0, cos^{-1}(|z'|/|r|), tan^{-1}(y'/x'))] \tag{4}
$$

 $r_0>R_{cut}/\pi$ → so that the 4D surface contains all the information from 3D spherical region inside the cutoff, include radial dimension

 **$\text{rotation angel}\ \omega \equiv \theta_0 \ \text{about some axis} \ n(\Theta \equiv \theta, \Phi \equiv \phi)$**  

**Condition:** 

Ref.[5] session 1.4.2, p.23

(1) $0 \leq \theta_0 \leq \pi$

(2) $0\leq \theta \leq \pi$

(3) $0 \leq \phi \leq 2\pi$.

 ************Note:  $\phi = tan^{-1}(\frac{y'}{x'})$, $x', y'$  has negative values → need to convert range of angle from $[-\pi, \pi] \ \text{to}\ [0, 2\pi]$**

(4) $R_{cut}$ depends on the chemical identities of both the neighbor atom and center atom (at a distance that is less force affected on both of them?)  

**************Compute $\theta_0, \ \theta, \ \phi$**

$$
\theta_{0} = \theta_0^{max}\frac{|r_{ik}|}{|r_0|}=\pi\frac{\sqrt{x_k^2+y^2_k+z_k^2}}{R_{cut}} \tag{5}
$$

where $\theta_0^{max}= \pi$ 

→ excluded points south  of latitude $\theta_0=r_{frac0}\pi$

```python
#theta_0
r_ik_array = df_ik['r_ik'].to_numpy() #r_ik from selected neighbors
r_0_array = np.full((r_ik_array.shape),R_cut)
theta_0_array = np.pi*(np.divide(r_ik_array,r_0_array))
print (r_0_array)
#theta
Z_k_abs_array = np.abs(df_ik['Z_k'].to_numpy())
theta_array = np.arccos(np.divide(Z_k_abs_array,r_ik_array))
#phi
X_k_array = df_ik['X_k'].to_numpy()
Y_k_array = df_ik['Y_k'].to_numpy()
phi_array = np.arctan(np.divide(Y_k_array, X_k_array))
#convert angle to positive value between [0,2pi]
phi_array_convert = np.mod(phi_array, 2*np.pi)
for angle in phi_array_convert:
    if (angle >=2*np.pi) and (angle < 0):
        raise ValueError('phi angle in between 0 and 2pi')
print(phi_array_convert)
print(phi_array)
#Replace NaN with 0: (code will have error for invalid value center atom values 0/0)_
df_ik['theta_0'] = theta_0_array
df_ik['theta_0'] = df_ik['theta_0'].replace(np.nan,0)
df_ik['theta'] = theta_array
df_ik['theta'] = df_ik['theta'].replace(np.nan,0)
df_ik['phi'] = phi_array_convert
df_ik['phi'] = df_ik['phi'].replace(np.nan,0)
print(df_ik.drop(['X', 'Y', 'Z', 'atom_type'], axis =1))
```

![Output 5. Data frame for possible neighbor atoms k  around a center atom i](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/Screen_Shot_2022-10-26_at_10.38.48_PM.png)

Output 5. Data frame for possible neighbor atoms k  around a center atom i

## 2.2 Compute **bispectrum component for input values $(j_1, j_2, j)$**

---

---

### Clebsch- Gordan coefficient

Ref.[3], Eq.(5) and Ref.[1], Eq.(5)

$$
C^{jm}_{{j_1}{m_1}{j_2}{m_2}} C^{jm'}_{{j_1}{m'_1}{j_2}{m'_2}} \equiv H^{jmm'}_{{{j_1}{m_1}{m'_1}} ,{{j_2}{m_2}{m'_2}}} \tag{0}  
$$

**Definition**

Ref.[5] p.235

A Clebsch-Gordan coefficients are vector addition coefficients. They play an important role in decomposition of reducible representations of rotational….Let angular momentum $j_1$ and $j_2$ with projections on $m_1$  and $m_2$ on the quantization axis. represents **the probability amplitude that $j_1$ and $j_2$ are coupled into a resultant angular momentum $j$ with projection $m$.**

$$
(j_1,m_1)\otimes(j_2,m_2) \to (jm) \tag{1}
$$

$$
C\equiv\braket{j_1m_1j_2m_2\vert jm} \tag{2}
$$

$$
\braket{j_1m_1j_2m_2\vert jm}=(-1)^{j_1-j_2+m}\sqrt{2j+1}\begin{pmatrix}j_1&j_2&j \\ m_1&m_2&{-m}\end{pmatrix} \tag{3}
$$

- AUTHORS: Jens Rasch (2009-03-24)- initial version
    
    ![Screen Shot 2022-08-03 at 7.09.53 PM.png](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/Screen_Shot_2022-08-03_at_7.09.53_PM.png)
    

---

### a. **Compute $C_{j_1m_1j_2m_2}^{jm}, \ C_{j_1m_1j_2m_2}^{jm'}$ via its relation to Wigner 3-j symbols**

Use function from sympy: already implied condition for $j_1,j_2, j, m_1, m_2, m$

```python
class sympy.physics.quantum.cg.CG(j1,m1,j2,m2,j3,m3)
cg = CG(j1,m1,j2,m2,j,m)
cg = cg.doit()
cg_p = CG(j1,m1p,j2,m2p,j,mp)
cg_p = cg_p.doit()
```

**Condition:**

In accordance with vector addition rules $j_1+j_2=j$, unless the triangular conditions (triangular in equalities) are fulfilled

(1)$|j_1-j_2|\ \le j \le j_1+j_2$ 

$j= 0, \frac{1}{2},1,...\infty$  and $m, m' = -j,-j+1, ...,j-1,j$

(2)$m_1+m_2=m \ \text{and} \ j_1 +j_2=j$

(3)$|m_1| \le j_1, \ |m_2| \le j_2,\ |m| \le j$ 

(4) $j_1,j_2,j$ not exceed a positive integer $J \ \text{is}$

(5) $j_1+j_2-j$ not half integer

(6) $m_1,m_2,m$ are integer or half-integer (positive or negative) numbers

(7) $j_1,j_2,j$ are integer or half integer non negative numbers 

(8) $j_1 +m_1, \ j_2 +m_2, \ j +m, \ j_1+j_2+j$  are integer non-negative numbers

**Parameter:**

$j_1, j_2$ : Angular momenta of states 1 and 2

$j,m\ :$  Total angular momentum of $(j_1+j_2)$

$m_1, \ m_2, \ m:$  Eigenvalues w.r.t to anglar momentum $j_1, \ j_2, \ j$

$m_1', \ m_2', \ m':$ Eigenvalues w.r.t $j_1, \ j_2, \ j$ along rotated axis

---

### b. **Compute coupling coefficient function $H^{jmm'}_{{{j_1}{m_1}{m'_1}} ,{{j_2}{m_2}{m'_2}}}$**

```python
import sympy as sp
from sympy.physics.quantum.cg import CG
def getCoeffH(j1,j2,j,m1,m2,m,m1p,m2p,mp):
    '''
    This function calculate coupling coefficient H via computing
    the Clebsch-Gordan coefficient for cg(j1,m1,j2,m2,j,m)
    and cg_p(j1,m1p,j2,m2p,j,mp)
    Parameters:
        j1: angular momentum 1
        j2: angular momentum 2
        j: total angular momentum (j1+j2)
        m1: eigenvalue of angular momentum j1
        m2: eigenvalue of angular momentum j2
        m: eigenvalue of angular momentum j
        m1p: eigenvalue of j1 along rotated axis
        m2p: eigenvalue of j2 along rotated axis
        mp: eigenvalue of j along rotated axis
    Returns:
        Coupling coefficient H(j1,j2,j,m1,m2,m.m1p,m2p,mp)
    ====================Reference=========================
    [1] Thompson, Swiler, Trott, Foiles, Tucker,
        Spectral neighbor analysis method for automated generation of quantum-accurate interatomic potentials (2015)
    [5] Chapter 8  D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
        Quantum Theory of Angular Momentum (1988)
    '''
    cg = CG(j1,m1,j2,m2,j,m)
    cg = cg.doit()
    cg_p = CG(j1,m1p,j2,m2p,j,mp)
    cg_p = cg_p.doit()
    H_coeff = (cg)*(cg_p)
		H = N(H_coeff)
    return H
```

**Example:**

In [ ]

```python
j1,j2,j,m1,m2,m,m1p,m2p,mp = (3/2,1/2,1,3/2,-1/2,1,3/2,1/2,2)
H = getCoeffH(j1,j2,j,m1,m2,m,m1p,m2p,mp)
print(H)
```

Out [ ]

```python
0
```

---

### c. **Compute expansion coefficients density function $u^j_{m,m'}, \ u^{j_1}_{m_1,m_2}, \ u^{j_2}_{m_1,m_2}$**

**Expansion coefficient for the partial density of neighbors of element $\mu$** 

**Definition**

$$
 u^j_{m,m'} = \braket{U^j_{m'm}\vert\rho}\tag{0}
$$

$$
\begin{align*} u_{j,m,m'}^{\mu} &= w^{self}_{\mu_i \mu} U_{j,m,m'}(0,0,0) \\ & + \displaystyle\sum_{r_{ik}<R_{cut}^{\mu_i\mu_k}} f_c(r_{ik};R_{cut}^{\mu_i\mu_k})w_{{\mu}_k}U_{m,m'}^{j}(\theta_{0},\theta, \phi) \delta_{\mu\mu_k} \end{align*} \tag{1}
$$

- $U_{m,m'}^{j}(\theta_{0},\theta, \phi)$
    
     Ref[5] Varshalovich, D A, Quantum Theory of Angular Momentum. 1988,  session 4.5 p.80
    
    $$
    U^j_{mm'}(\theta_0,\theta,\phi)=\displaystyle\sum_{m''=-j} ^{j} D^j_{mm''}(\phi,\theta,-\phi) e^{-im''\theta_0}D_{m''m'}(\phi, -\theta,- \phi) \tag{4}
    $$
    
    ![Screen Shot 2022-09-15 at 3.24.55 PM.png](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/Screen_Shot_2022-09-15_at_3.24.55_PM.png)
    
    **Condition: $\text{rotation angel}\ \omega \equiv \theta_0 \text{about some axis} \ n(\Theta \equiv \theta, \Phi \equiv \phi)$**  
    
    (1) $0 \leq \theta_0 \leq \pi$
    
    (2) $0\leq \theta \leq \pi$
    
    (3) $0 \leq \phi \leq 2\pi$  
    
    **Explicit forms of the Wigner D-functions**
    
    Ref.[5] Varshalovich Eq. (1), Session 4.3. page 76
    
    $D^j_{mm'}(\alpha,\beta, \gamma)$  represents as a product of three functions each of which depends on one arguement $\alpha, \beta$ or $\gamma$,
    
    $$
    D^j_{mm'}(\alpha, \beta,\gamma)= \ e ^{-im\alpha}d^j_{mm'}(\beta)e^{-im'\gamma} \tag{5}
    $$
    
    which we can evaluate the Wigner D matrix elements of a rotation $D^j_{mm'}(\phi,\theta,-\theta)$  and  using sympy package
    
    $$
     D^j_{mm''}(\phi,\theta,-\phi) e^{-im''\theta_0}D_{m''m'}(\phi, -\theta,- \phi) \tag{6}
    $$
    
    **Compute $U_{m,m'}^{j}(\theta_{0},\theta, \phi)$**
    
    ```python
    import sympy as sp
    from sympy.physics.quantum.cg import CG
    from sympy.physics.wigner import wigner_d
    from sympy.physics.quantum.spin import Rotation
    from sympy import *
    ```
    
    ```python
    def getRotationalMatrixU(j,m,mp,mpp,theta_0,theta,phi):
        '''
        Parameters:
            j: integer/half integer number
               Total angular momentum
            m: integer/half integer number
               Eigenvalue of angular momentum along rotated axis
            mp:
            mpp:
            theta_0: fist angle of rotation [0,pi]
            theta: second angle of rotation [0,pi]
            phi: third angle of rotation [0,2pi]
            rotational od a coordinate system through an angle theta_0
            about an axis n(theta,phi)
        Returns: Rotational matrix U
        ==========================Reference==================================
        [5] Chapter 4  D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
            Quantum Theory of Angular Momentum (1988)
        '''
        mpp_list = np.linspace(-j, j, int(2 * j + 1))
        U = 0
        for mpp in mpp_list:
            rot1 = Rotation.D(j, m, mpp, phi, theta, - phi)
            rot2 = Rotation.D(j, mpp, mp, phi, -theta, -phi)
            Dj_m_mpp = rot1.doit()
            Dj_mpp_mp = rot2.doit()
            Um_mp = Dj_m_mpp * (exp(-I * mpp * theta_0)) * Dj_mpp_mp
            Um_mp = nsimplify(Um_mp)
            U += N(Um_mp)
        return U
    ```
    
    **Example**
    
    → because $U^j_{mm'}(\theta_0,\theta,\phi)$ define as eigenfunctions Ref.[3], session 4.5.3, page 82. It makes sense that the output is complex numbers which represents a rotation 
    
- $f_c(r_{ik};R_{cut})$
    - $f_c$ is the cutoff function to ensure smooth decay for the neighbor atomic density to zero cutoff radius $R_{cut}^{\mu_i\mu_k}$ Ref[2]
        
        $$
        f_c(r) = \frac{1}{2}(cos(\pi\frac{r_{ik}-r_{min0}}{R^{\mu_i\mu_k}_{cut}-r_{min0}})+1),r_{ik}\leq R_{cut}^{\mu_i\mu_k}\\=0, r_{ik}>R^{\mu_i\mu_k}_{cut} \hspace{2.6cm} \tag{2}
        $$
        
    
    Ref.[7]
    
    $$
    f_c(r) = \frac{1}{2}(cos(\pi\frac{r_{ik}-r_{min0}}{R_{ik}-r_{min0}})+1),\leq R_{ik}
    $$
    
    - $r_{ik}$: distance from neighbor atom to center atom.
    - $R_i, R_k$ (center, neighbor): atom radius, measures from center of nucleus to the outermost isolated electron.
    - $r_{min0}$: parameter in distance to angle conversion (distance units), choose $0$
        
        **Compute $f_c(r_{ik},R_{cut}^{ik})$ cut off function**
        
        ```python
        def getCutoffFunction (r_ik,r_min0, R_cut):
            '''
            Parameter:
                r_ik: 1D array dimension [k,], distance from center atom  to k neighbor atoms
                r_min0: parameter in distance to angle conversion(distance unit)
                R_cut: cut off radius
            Returns:
                array, dimension [k,]: cut of function f_cut(r_ik) where r_ik < R_cut
            ==========================References==================================
            [3] A. Bartok, M. C. Payne, K. Risi, G. Csanyi, Gaussian approximation
                potentials: the accuracy of quantum mechanics, without the electrons (2010)
            [7] LAMMPS, Compute sna/atom command,
                https://docs.lammps.org/compute_sna_atom.html
            '''
            R_cut_array = np.full((r_ik_array.shape), R_cut)
            f_cut = (1 / 2) * (np.cos(np.pi * (np.divide(r_ik_array - r_min0, R_cut_array - r_min0))) + 1)
            return f_cut
        ```
        
- $\mu$ is chemical element
- $w_{\mu_i\mu}^{self}$ is self contribution ( either 1 or 0) but we can ignore this since $U^j_{mm'}(0,0,0) = 0$
- $\omega_{\mu_k}$ are the coefficients that are dimensionless weights that are chosen to distinguish atoms of different types, while the central atom is arbitrarily assigned a unit weight. Ref.[1]
- $\delta_{\mu \mu_k}$ is the Dirac delta function, indicates only neighbor atom of element $\mu$ contribute to partial density $\rho^{\mu}$ by requiring that the partial densities sum to the total density used in Req.(7)
    
    $$
    \rho(r)=\displaystyle\sum_{\mu =1}^{N_{elem}}\rho^{\mu}(r) \tag{3}
    $$
    
    ********************************************Compute $u^j_{m,m'}$ function**
    
    - **************Example**************
        
        In [ ]
        
        ```python
        #INPUT VALUES
        j,m,mp = 3,2,3
        #array for weight coefficient w.r.t to atome type
        
        #Cutoff function
        def getCutoffFunction(r_ik_array, r_min0, R_cut):
            R_cut_array = np.full((r_ik_array.shape),R_cut)
            f_cut = (1/2)*(np.cos(np.pi*(np.divide(r_ik_array-r_min0, R_cut_array-r_min0)))+1)
            return f_cut
        f_cut_arr = getCutoffFunction(r_ik_array,0, R_cut)
        print (f_cut_arr)
        #U_immp
        from methods import *
        U_jmmp=[]
        for theta_0, theta, phi in zip(theta_0_array, theta_array, phi_array):
            U = getRotationalMatrixU(j, m,mp,theta_0, theta, phi)
            U_jmmp.append(U)
        U_jmmp_array = np.array(U_jmmp, dtype='complex')
        U_jmmp_arr = np.where(np.isnan(U_jmmp_array), 0, U_jmmp_array)
        print(U_jmmp_arr)
        #Create array dimension [k,] for w_ik,
        #because all atoms are Si so we set weight coeffient =1 if not w_ik=w_i/w_k
        w_ik_arr = np.full((r_ik_array.shape),1)
        #delta function delta=1 if i and k has the same element type, if not delta =0
        delta = np.full((r_ik_array.shape),0)
        delta_arr = np.where(df_ik['atom_type']==df_ik['atom_type'].iloc[0],1,delta)
        print(delta_arr)
        u_jmmp = np.dot((f_cut_arr*U_jmmp_arr),(w_ik_arr*delta_arr))
        print(u_jmmp)
        ```
        
        Out [ ]
        
        ```python
        (-0.43441760999036916+0.3597729040449244j)
        ```
        
    
    ```python
    def getDensityFunction_u(j,m,mp,w_ik_array, delta_array,r_ik_array, r_min0, R_cut, theta_0_array, theta_array, phi_array):
        '''
        Args:
            j: angular momentum
            m: eigenvalue of angular momentum
            mp: eigenvalue of j along rotated axis
            w_ik_arr: array for the coefficients that are dimensionless weights that are chosen to distinguish atoms
                        of different types, while the central atom is arbitrarily assigned a unit weight, dimensin (1,k)
            delta_arr: array for the Dirac delta function, indicates only neighbor atom of element the same as center atom
                        contribute to partial density,  dimension (1,k)
            r_ik_array: array for distance from center atom to neighbor atom, dimension (1,k+1), k is number of neighbor atoms
                        in cutoff radius, array include center atom as well
            r_min0: number, parameter in distance to angle conversion (distance units), choose
            R_cut: number, cutoff radius
            theta_0_array: array for theta_0 angel (fist angle of rotation [0,pi])
                        of neighbor atoms in reference frame of center atom, dimension (1, k+1)
            theta_array: array for theta angel ( second angle of rotation [0,pi])
                        of neighbor atoms in reference frame of center atom, dimension (1, k+1)
            phi_array: array for phi angel (third angle of rotation [0,2pi])
                        of neighbor atoms in reference frame of center atom, dimension (1, k+1)
        Returns: expansion coefficients density function u_jm_mp
        '''
        R_cut_array = np.full((r_ik_array.shape), R_cut)
        f_cut_arr = (1 / 2) * (np.cos(np.pi * (np.divide(r_ik_array - r_min0, R_cut_array - r_min0))) + 1)
        U_jmmp = []
        for theta_0, theta, phi in zip(theta_0_array, theta_array, phi_array):
            mpp_list = np.linspace(-j, j, int(2 * j + 1))
            U = 0
            for mpp in mpp_list:
                rot1 = Rotation.D(j, m, mpp, phi, theta, - phi)
                rot2 = Rotation.D(j, mpp, mp, phi, -theta, -phi)
                Dj_m_mpp = rot1.doit()
                Dj_mpp_mp = rot2.doit()
                Um_mp = Dj_m_mpp * (exp(-I * mpp * theta_0)) * Dj_mpp_mp
                Um_mp = nsimplify(Um_mp)
                U += N(Um_mp)
            U_jmmp.append(U)
        U_jmmp_array = np.array(U_jmmp, dtype='complex')
        U_jmmp_arr = np.where(np.isnan(U_jmmp_array), 0, U_jmmp_array)
        u_jmmp = np.dot((f_cut_arr*U_jmmp_arr),(w_ik_arr*delta_arr))
        return u_jmmp
    ```
    
    **Comment:** 
    
    Because  $U_{j,m,m'}(0,0,0) = 0$  so we ignore self contribution $w_{\mu_i\mu}^{self}$. 
    
    We simplify: 
    
    $$
    \begin{align*} u_{j,m,m'}^{\mu} &= w^{self}_{\mu_i \mu} U_{j,m,m'}(0,0,0) \\ & + \displaystyle\sum_{r_{ik}<R_{cut}^{\mu_i\mu_k}} f_c(r_{ik};R_{cut}^{\mu_i\mu_k})w_{{\mu}_k}U_{m,m'}^{j}(\theta_{0},\theta, \phi) \delta_{\mu\mu_k} \end{align*} \tag{1}
    $$
    
    to
    
    $$
    u_{j,m,m'}^{\mu} =  \displaystyle\sum_{r_{ik}<R_{cut}^{\mu_i\mu_k}} f_c(r_{ik};R_{cut}^{\mu_i\mu_k})w_{{\mu}_{ik}}U_{m,m'}^{j}(\theta_{0},\theta, \phi) \delta_{\mu\mu_k} 
    $$
    
    ![Output 6. Visualize input values as columns for calculating the summation of the expansion coefficient density function u (j=3, m=2, m’=3) ](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/Screen_Shot_2022-10-31_at_3.49.09_PM.png)
    
    Output 6. Visualize input values as columns for calculating the summation of the expansion coefficient density function u (j=3, m=2, m’=3) 
    

### Compute bispectrum components

**********************Example input : $j_1, \ j_2, \ j = 1, \ 2, \ 3$**

$$
B_(j_1,j_2,j) = \displaystyle\sum_{{m,m'= -j}}^j  \left(u^{\smash{j}}_{m,m'}\right)^*\displaystyle\sum_{{m_1,m'_1 =-j_1}}^{j_1} \displaystyle\sum_{{m_2,m'_2 = -j_2}}^{j_2} H^{jmm'}_{{{j_1}{m_1}{m'_1}} ,{{j_2}{m_2}{m'_2}}}  u^{j_1}_{{m_1},{m_1'}} u^{j_2}_{{m_2},{m_2'}} 
$$

- $C^{jm}_{{j_1}{m_1}{j_2}{m_2}} C^{jm'}_{{j_1}{m'_1}{j_2}{m'_2}} \equiv H^{jmm'}_{{{j_1}{m_1}{m'_1}} ,{{j_2}{m_2}{m'_2}}}$
- $\displaystyle\sum_{{m,m'= -j}}^j  \left(u^{\smash{j}}_{m,m'}\right)^*$
- $\displaystyle\sum_{{m_1,m'_1 =-j_1}}^{j_1} \displaystyle\sum_{{m_2,m'_2 = -j_2}}^{j_2} H^{jmm'}_{{{j_1}{m_1}{m'_1}} ,{{j_2}{m_2}{m'_2}}}  u^{j_1}_{{m_1},{m_1'}} u^{j_2}_{{m_2},{m_2'}}$

****************************Question:** 

- Does $(j_1, \ j_2, \ j)$ depend on the center atom type?
- Does that means  we need to compute possible values for $m_1, \ m_2, \ m$, and $m_1', \ m_2', \ m'$? which needs to satisfy the condition listed in part 2.2.a?
    - Doesn’t have to, when using sympy to calculate coupling coefficient $H$, invalid input values $(j,j_1,j_2,m,m_1,m_2,m',m_1', m_2')$ return to 0 so we would just skip set that $H=0$

[https://stackoverflow.com/questions/533905/get-the-cartesian-product-of-a-series-of-lists](https://stackoverflow.com/questions/533905/get-the-cartesian-product-of-a-series-of-lists)

**Idea:** create possible set $(j_1,j_2,j,m_1,m_2,m,m'_1,m_2',m')$ from input $(j_1, j_2,j)$

```python
#create input set [j1,j2,j,m1,m2,m,m1p,m2p,mp]
j = 3
j1 = 1
j2 = 2
m = np.linspace(-j, j, int(2 * j + 1)).tolist()
mp = np.linspace(-j, j, int(2 * j + 1)).tolist()
m1 = np.linspace(-j1, j1, int(2 * j1 + 1)).tolist()
m1p = np.linspace(-j1, j1, int(2 * j1 + 1)).tolist()
m2 = np.linspace(-j2, j2, int(2 * j2 + 1)).tolist()
m2p = np.linspace(-j2, j2, int(2 * j2 + 1)).tolist()
```

![IMG_2050.jpg](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/IMG_2050.jpg)

```python
B=0
for i in itertools.product(m1,m2,m,m1p,m2p,mp):
    m1, m2, m, m1p, m2p, mp = i 
    j,j1,j2=3,1,3
    H = getCoeffH(j1,j2,j,m1,m2,m,m1p,m2p,mp)
    if H==0:
        pass
    else:
        u_jmmp = getDensityFunction_u(j, m, mp, w_ik_arr, delta_arr, r_ik_array, 0,
                                      R_cut, theta_0_array, theta_array,phi_array)
        u1_j1m1m1p = getDensityFunction_u(j1, m1, m1p, w_ik_arr, delta_arr, r_ik_array, 0,
                                          R_cut, theta_0_array,theta_array, phi_array)
        u2_j2m2m2p = getDensityFunction_u(j2, m2, m2p, w_ik_arr, delta_arr, r_ik_array, 0,
                                          R_cut, theta_0_array,theta_array, phi_array)
        B_each = np.conj(u_jmmp) * H * (u1_j1m1m1p) * (u2_j2m2m2p)
        B = N(B_each)
				B_sum += B
print (B_sum)
```

Out [ ]

```python
3.89033918818112 + 4.16333634234434e-17*I
```

Con: running time too long

---

# References

---

[1] `*Thompson, Swiler, Trott, Foiles, Tucker,*` ***Spectral neighbor analysis method for automated generation of quantum-accurate interatomic potentials***  (2015) [https://www.osti.gov/biblio/1426894](https://www.osti.gov/biblio/1426894) 

[2]  `*J. K. Mason,` **The relationship of the hyperspherical harmonics to SO(3), SO(4) and orientation distribution functions***, Acta Cryst A65, 259 (2009)

[https://libraryh3lp.com/file/8j46b06z956ega%40web.libraryh3lp.com/1665598079.pdf?t=6sfhoB1XgCcztCUZBmTEDs](https://libraryh3lp.com/file/8j46b06z956ega%40web.libraryh3lp.com/1665598079.pdf?t=6sfhoB1XgCcztCUZBmTEDs)

[3] *`A. Bartok, M. C. Payne, K. Risi, G. Csanyi,`**Gaussian approximation potentials: The accuracy of quantum mechanics, without the electrons*** (2010) [https://arxiv.org/pdf/0910.1019.pdf](https://arxiv.org/pdf/0910.1019.pdf)

[4] `*A V Meremianin,*`  ***Multipole expansions in four-dimensional hyperspherical harmonics,*** J. Phys. A: Math. Gen. 39 3099 (2006) [https://iopscience.iop.org/article/10.1088/0305-4470/39/12/017/pdf?casa_token=YfzEUY2g4jwAAAAA:bqMuUwTRpQXDQEpCSvvwmlFYX6hi0xis-vCqLThjemvfDObHjP7XZw28oexUMra9FGg7AV1FVKvhtzZJn28g](https://iopscience.iop.org/article/10.1088/0305-4470/39/12/017/pdf?casa_token=YfzEUY2g4jwAAAAA:bqMuUwTRpQXDQEpCSvvwmlFYX6hi0xis-vCqLThjemvfDObHjP7XZw28oexUMra9FGg7AV1FVKvhtzZJn28g)

[5] *`D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,` **Quantum Theory of Angular Momentum*** (1988) [https://library.oapen.org/handle/20.500.12657/50493](https://library.oapen.org/handle/20.500.12657/50493)

[6]  *`M. A. Cusentino,  M. A. Wood, A. P. Thompson,`* ***Explicit Multi-element Extension of the Spectral Neighbor Analysis Potential for Chemically Complex Systems,*** J. Phys. Chem. A,**[](https://pubs.acs.org/action/showCitFormats?doi=10.1021%2Facs.jpca.0c02450&href=/doi/10.1021%2Facs.jpca.0c02450)** 124, 26, 5456–5464 (2020) [https://doi.org/10.1021/acs.jpca.0c02450](https://doi.org/10.1021/acs.jpca.0c02450)

[7] `*LAMMPS,*` ***Compute sna/atom command,*** [https://docs.lammps.org/compute_sna_atom.html](https://docs.lammps.org/compute_sna_atom.html)

*[8] `Kyushin, S., Kurosaki, Y., Otsuka, K. et al.` ****Silicon–silicon π single bond**** Nat Commun* **11**, 4009 (2020). https://doi.org/10.1038/s41467-020-17815-z

- extra
    
    [8] Complex strengthening mechanisms in the NbMoTaW multi-principal element alloy [https://doi.org/10.1038/s41524-020-0339-0](https://doi.org/10.1038/s41524-020-0339-0)
    
    [9] Multiple expansion in four-dimensional hyperspherical harmonic [https://arxiv.org/pdf/math-ph/0510080.pdf](https://arxiv.org/pdf/math-ph/0510080.pdf)
    
    [10] Performance and cost assessment of machine learning interatomic potentials
    
    chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/[https://materialsvirtuallab.org/pubs/10.1021_acs.jpca.9b08723.pdf](https://materialsvirtuallab.org/pubs/10.1021_acs.jpca.9b08723.pdf)
    
    [11] LAMMPS Feature and Capabilities [https://www.osti.gov/servlets/purl/1315254](https://www.osti.gov/servlets/purl/1315254)
    
    [12] LAMMPS website [https://www.lammps.org/cite.html](https://www.lammps.org/cite.html)
    
    [13] Journal of Mathematical Physics: A calculation of SU(4) Clebsch- Gordan
    
    [https://umkc-illiad-oclc-org.proxy.library.umkc.edu/illiad/illiad.dll](https://umkc-illiad-oclc-org.proxy.library.umkc.edu/illiad/illiad.dll)
    
    [14] Determination of Clebsch- Gordan coefficients by matrix diagonalization [https://aip.scitation.org/doi/pdf/10.1063/1.168419](https://aip.scitation.org/doi/pdf/10.1063/1.168419)
    
    [15] Efficient storage scheme for precalculated wigner 3j, 6j and gaunt coeffcients
    
    [https://www.theochem.ru.nl/files/local/sjsc-25-1416-2004.pdf](https://www.theochem.ru.nl/files/local/sjsc-25-1416-2004.pdf)
    
     [16] Theory of Angular Momentum
    
    [https://kgut.ac.ir/useruploads/1505647831850hcd.pdf](https://kgut.ac.ir/useruploads/1505647831850hcd.pdf)
    
    [17] Gaussian Approximation Potential an interatomic potential derived from first principles Quantum Mechanics [https://arxiv.org/pdf/1003.2817.pdf](https://arxiv.org/pdf/1003.2817.pdf)
    

**GitHub**

[1][https://github.com/zichunhao/WignerD/blob/main/cg/clebsch_gordan_coefficients.py](https://github.com/zichunhao/WignerD/blob/main/cg/clebsch_gordan_coefficients.py)

[2] LAMMPS [https://github.com/cesmix-mit/LAMMPS.jl/blob/main/examples/snap.jl](https://github.com/cesmix-mit/LAMMPS.jl/blob/main/examples/snap.jl)

[3] Wigner-D [https://github.com/sympy/sympy/blob/c79d74dafb21d631f54ba82d5ddecc56dc9efaa3/sympy/physics/quantum/spin.py#L699-L904](https://github.com/sympy/sympy/blob/c79d74dafb21d631f54ba82d5ddecc56dc9efaa3/sympy/physics/quantum/spin.py#L699-L904)

## Contact

---

![download.png](Compute%20Bispectrum%20Components%20fe15e8c590c94d438538084914489ecb/download.png)

### Duong Hoang

Email address: dtht2d@umsystem.edu

****UMKC Computational Physics Group****

**School of Science and Engineering**

**Division of Energy, Matter, and Systems**
