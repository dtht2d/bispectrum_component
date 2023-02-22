# **Compute Bispectrum Components**
---
## Motivation

This research aims to leverage Python scripts to compute bispectrum components, which will then be used as descriptors to train a machine learning model that is capable of predicting the electronic structure of materials. A comprehensive grasp of the electronic structure of materials is fundamental for gaining insight into their behavior at the molecular and nanoscale level. The arrangement of electrons in atoms or molecules is what determines their physical and chemical properties, and the development of predictive models for the electronic structure of materials can facilitate the optimization of their properties.

This study aims to further understand the correlation between atomic structure and its local environment, and to develop models that can explore the electronic properties of novel materials and enhance the characteristics of existing materials. Electronic structure is a critical factor in several properties of materials such as how atoms bond to form molecules, how they interact with light, and the material's electrical conductivity, thermal conductivity, and mechanical properties. A comprehensive comprehension of electronic structure is necessary for the creation of new materials with superior properties and the optimization of the performance of existing materials.
Ultimately, the development of precise electronic structure prediction models can facilitate the design and engineering of advanced materials for various applications, including electronics, energy storage, catalysis, and drug discovery. Thus, the progress of our knowledge of electronic structure is essential for the development of materials science and the invention of new technologies that can benefit society.

---
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
