# **Compute Bispectrum Components**
---
## Objective

The aim of this research is to employ Python scripts to compute bispectrum components, which will then be utilized as descriptors for training a machine learning model to predict material force fields. This study aims to further our comprehension of the relationship between atomic structure and force fields, and to apply this knowledge towards the development of predictive models for material behavior at the molecular and nanoscale level. Material force fields are mathematical models that simulate the interactions between atoms and molecules in a material and provide valuable information on its structure, thermodynamics, and mechanical properties. This information is critical for various applications, including the design and optimization of new materials, the prediction of material behavior in varying environments, and the advancement of manufacturing processes. Through the study of material force fields, scientists and engineers can gain deeper insight into the underlying physical principles governing material behavior, leading to the improvement of existing materials and the creation of new materials with superior properties for diverse applications.

---
## Outline
**Computational steps**

**1. Data prep and parameter**

**a. Data prep**

- Read CIF file as data frame

**b. Define position of neighbor atoms $k$ relative to a central atom $i$  is a point within the 3D ball of radius**

- Estimate list of potentially atoms in the center cell
- Choose a center atom $i$  from the list and get its coordinate  ( $x_i,y_i,z_i$ )
- Re-calculate coordinate for all atoms in the cell w.r.t new reference of frame  (origin at ( $x_i, y_i, z_i$ ) )
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
 
