# Abstract 
First-principles electronic structure calculations based on density functional theory (DFT) are well-known to have a high computational cost that scales algorithmically as $O(N^3)$, where $N$ is the number of electrons. Reducing that cost is a key goal of the computational materials physics community and machine learning (ML) is viewed as a potential tool for that task. However, ML model training requires carefully selected input descriptors for training. This resporsitory presents a computer program codes (bispectrum/methods/calc) that is designed to automate the generation of local atomic environment descriptors for single element systems. The descriptors may be used for training neural networks to predict the set of electronic potential function coefficients, $\{A_i\}$, that are used within the DFT based orthogonalized linear combination of atomic orbitals (OLCAO) method. Predicting the potential function coefficients rather than explicitly computing them will reduce the computational cost of a typical OLCAO calculation by at least one order of magnitude and possibly two. We explore additional research directions by connecting the gap between descriptors used in computer vision and those employed in electronic structure predictions. This approach opens up possibilities for cross-disciplinary knowledge and techniques from computer vision, which may further improve the accuracy and efficiency of electronic structure calculations.
#Premilinary Results
![Example](https://github.com/dtht2d/bispectrum_component/blob/main/plots/various-si-models-testing/cif-files/I4-si.png)

![plot](https://github.com/dtht2d/bispectrum_component/blob/main/figures/I4-si-origin.png)

