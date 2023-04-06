# **Compute Bispectrum Components**
## Introduction
---

In the field of materials science, machine learning and computational physics have revolutionized the process of discovering and developing new materials. Accurately predicting electronic structure is crucial in designing materials with desired properties. Although Density Functional Theory (DFT) is a popular computational model, it faces challenges when dealing with transition metals, rare earth elements, and strong correlation effects. Therefore, current research aims to improve DFT's accuracy and applicability by developing new functionals and techniques.

Objective: This research aims to develop a new machine learning framework and additional descriptors such as bispectrum components, Coulomb matrix eigenvalues (CMEs), and Atom-centered Symmetry Functions (ACSFs) to enhance the accuracy and applicability of predicting electronic structure. The focus of this repository is to optimize the code for computing bispectrum components using Python and use them as input features for a machine learning model.

Significance: Predicting electronic structure enhances our understanding of the correlation between atomic structure and local environment, which can facilitate the development of advanced materials and new technologies. Bispectrum components are a promising descriptor for capturing essential features of electronic structure that traditional spectral descriptors may miss.

Measure of Success: The success of this project will be the development of Python scripts that can automatically generate bispectrum coefficients from crystal structure datasets. The accuracy of the calculations will be tested by comparing predicted electronic structure results to experimental or theoretical results. This research has the potential to contribute to the field of materials science by providing more accurate and efficient methods for predicting electronic structure and designing new materials. Additionally, I am interested in utilizing their research to develop new mediums and techniques for creative expression.

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
---
**3. Code Optimization**

## Future work:
Propose new Machine Learning framework:

  - What type of Neural Networks would work?
  - And why?
  
Compute other descriptors for input training: 

  - Atom-centered Symmetry Functions (ACSFs): descriptor of the local environment of each atom, these models can capture the complex relationships between the atomic positions and the electronic structure of the molecule, leading to accurate predictions of electronic properties.
  - Coulomb matrix eigenvalues (CMEs) 
  - Bag of bond 

