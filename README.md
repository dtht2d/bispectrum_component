# **Compute Bispectrum Components**
---
## Abstract 
The development of a Python script to calculate bispectrum components and integrate them into the Orthogonalized Linear Combination of Atomic Orbitals (OLCAO) package is the focus of my research. Bispectrum components have proven to be a valuable and promising tool in predicting material properties, as they can capture subtle variations in the local atomic environment that are difficult to predict using other methods like neural networks and graph neural networks. By incorporating bispectrum components into machine learning frameworks, my research aims to more accurately predict material properties and electronic structures influenced by nearby atoms.

To evaluate the significance of this research, we have optimized the running time of the bispectrum calculation and developed new custom functions that don't rely on third-party libraries. The current work has calculated the bispectrum component, but the running time took a long time. Further work will involve code optimization of the bispectrum class function, visualizing the relationship between bispectrum components and the band limit ${(j,j_1,j_2)}$, and studying the expansion density function in hyper-spherical harmonic.

One challenge in this research is defining a suitable cut-off radius for the cut-off function to avoid neglecting the interaction between the center atom and its defined neighbor. Additionally, understanding the relationship between quantum numbers and the density coefficient function is crucial for accurately describing particle behavior in higher-dimensional space. This lack of understanding can limit the development of methods to predict the material's behavior based on the bispectrum, underscoring the urgent need for further research in this area to realize the full potential of the hyper-spherical harmonic expansion of density function.

In future work, we will compare bispectrum components to other descriptors for input training, including Atom-centered Symmetry Functions (ACSFs), Coulomb matrix eigenvalues (CMEs), and Bag of Bond. These models can capture complex relationships between the atomic positions and the electronic structure of the molecule, leading to accurate predictions of electronic properties. Overall, this research contributes to the ongoing effort to develop new and improved machine learning frameworks for predicting electronic structure materials with desirable properties and advancing the field of material science research.

Objective
---

**Stage 1** : Optimize function to calculate bispectrum components from given $j_1, j_2, j$

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

**3. Code Optimization**

---
## Future work:
    - Visualize relationship between bispectum component and band limit (j,j1,j2) 

    - Study expansion density function in hyper-spherical harmonic

    - Compare bispectrum to other descriptors for input training

        - Atom-centered Symmetry Functions (ACSFs): descriptor of the local environment of each atom, these models can capture the complex relationships between the atomic positions and the electronic structure of the molecule, leading to accurate predictions of electronic properties.​

        - Coulomb matrix eigenvalues (CMEs) 

        - Bag of bond 

