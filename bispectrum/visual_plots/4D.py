import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Generate points on the 4D unit sphere
num_points = 1000
theta = np.linspace(0, 2 * np.pi, num_points)
phi = np.linspace(0, np.pi, num_points)
u = np.linspace(0, 2 * np.pi, num_points)
v = np.linspace(0, np.pi, num_points)

x = np.outer(np.sin(theta), np.sin(phi)) * np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(theta), np.sin(phi)) * np.outer(np.sin(u), np.sin(v))
z = np.outer(np.cos(theta), np.sin(phi)) * np.outer(np.ones(num_points), np.cos(v))
w = np.outer(np.cos(theta), np.sin(phi)) * np.outer(np.ones(num_points), np.cos(v))

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the 4D sphere projection
ax.scatter(x, y, z, c=w, cmap='Spectral')

# Set plot properties
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('4D Unit Sphere Projection')

# Show the plot
plt.show()