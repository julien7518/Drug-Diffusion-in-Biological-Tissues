import numpy as np
import matplotlib.pyplot as plt

# Domain and grid
L = 1.0
N = 100
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)
X, Y = np.meshgrid(x, y)

D = 0.01  # diffusion coefficient
t = 0.5  # time
x0, y0 = 0.5, 0.5  # initial point source location

U = 1 / (4 * np.pi * D * t) * np.exp(-((X - x0)**2 + (Y - y0)**2) / (4 * D * t))

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, U, cmap='viridis', edgecolor='none')
ax.set_title(f'Analytical 2D Diffusion at t = {t} s (D = {D} cm²/s)')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_zlabel('u (mass per cm$^2$)')
plt.tight_layout()

plt.savefig("resources/figures/analytical_solution.png", dpi=300)
plt.show()