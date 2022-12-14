# **Compute Bispectrum Components**
---
## Outline
**Computational steps**

**1. Data prep and parameter**

**a. Data prep**

- Read CIF file as data frame

**b. Define position of neighbor atoms $k$ relative to a central atom $i$  is a point within the 3D ball of radius**

- Estimate list of potentially atoms in the center cell
- Choose a center atom $i$  from the list and get its coordinate $(x_i,y_i,z_i)$
- Re-calculate coordinate for all atoms in the cell w.r.t new reference of frame  (origin at $(x_i,y_i,z_i)$)
- Compute neighbors list in a chosen cutoff radius $R_{cut}$ with respect to new reference frame where its origin is at center atom $i$  location

**c. Map possible neighbor atoms on to the set of points $(\theta_0, \theta,\phi)$**

Compute $\theta_0, \ \theta, \ \phi$

**2. Compute bispectrum component for input values $(j_1, j_2, j)$**
1. Compute Clebsch- Gordan coefficient
$$ C^{jm}_{{j_1}{m_1}{j_2}{m_2}}, C^{jm'}_{{j_1}{m'_1}{j_2}{m'_2}} $$

2. Compute coupling coefficient
$$H^{jmm'}_{{{j_1}{m_1}{m'_1}} ,{{j_2}{m_2}{m'_2}}}$$

3. Compute  density coefficient
$$u^j_{mm'}, \mu_{m_1m_1'}^{j_1}, \mu_{m_2,m_2'}^{j_2}$$
---
 
