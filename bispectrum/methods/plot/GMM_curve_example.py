'''
This script plots the Gaussian curves used in the GMM method.
'''
import numpy as np
import matplotlib.pyplot as plt

# Generate x values
x = np.linspace(-4, 4, 1000)

# Compute the standard Gaussian distribution values
y1 = 1 / np.sqrt(2 * np.pi) * np.exp(-0.5 * x**2)

# Compute additional Gaussian distributions with different parameters
mu_values = [-2, 0, 1, 2, 4]
sigma_values = [0.5, 1, 1.5, 2, 2.5]
colors = ['red', 'green', 'blue', 'purple', 'orange']

plt.figure()
# Plot the curves
for i in range(5):
    y = 1 / (sigma_values[i] * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mu_values[i]) / sigma_values[i])**2)
    plt.plot(x, y, color=colors[i], label=f'mu={mu_values[i]}, sigma={sigma_values[i]}')

plt.title('Gaussian Distributions')
plt.xlabel('Distance')
plt.ylabel('Probability Density')
plt.savefig("/Users/duonghoang/Documents/GitHub/bispectrum_component/plots/GMM_example", dpi=300, bbox_inches='tight')
# Display the plot
plt.show()
