import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

L = 1.0              # domain size (cm)
N = 100              # number of grid points
dx = L / N           # spatial step
dy = dx              # assume dx = dy
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)
X, Y = np.meshgrid(x, y)

sigma = 0.05
T = 1.0              # total simulation time (s)

def analytical_solution(x, y, t, D):
    return 1 / (4 * np.pi * D * t) * np.exp(-((x - 0.5)**2 + (y - 0.5)**2) / (4 * D * t))

def run_simulation(D_x, D_y, lambda_decay, label):
    dt = 0.25 * min(dx**2 / D_x, dy**2 / D_y)
    Nt = int(T / dt)
    u = np.exp(-((X - 0.5)**2 + (Y - 0.5)**2) / (2 * sigma**2))
    results = [u.copy()]

    for n in range(Nt):
        u_new = u.copy()
        u_new[1:-1, 1:-1] = (
            u[1:-1, 1:-1] +
            dt * (
                D_x * (u[2:, 1:-1] - 2*u[1:-1, 1:-1] + u[:-2, 1:-1]) / dx**2 +
                D_y * (u[1:-1, 2:] - 2*u[1:-1, 1:-1] + u[1:-1, :-2]) / dy**2 -
                lambda_decay * u[1:-1, 1:-1]
            )
        )
        u_new[0, :] = u_new[1, :]
        u_new[-1, :] = u_new[-2, :]
        u_new[:, 0] = u_new[:, 1]
        u_new[:, -1] = u_new[:, -2]
        u = u_new
        if n % (Nt // 10) == 0:
            results.append(u.copy())

    if D_x == D_y and lambda_decay == 0:
        U_exact = analytical_solution(X, Y, T, D_x)
        rmse = np.sqrt(np.mean((U_exact - u)**2))
        print(f"[{label}] RMSE at t={T}s: {rmse:.6f}")
    else:
        print(f"[{label}] Convergence test skipped.")

    output_dir = f"resources/{label}"
    photo_dir = os.path.join(output_dir, "photos")
    os.makedirs(photo_dir, exist_ok=True)

    for i, u_plot in enumerate(results):
        plt.figure(figsize=(5, 4))
        plt.title(f"{label} - t = {round(i * dt * (Nt // 10), 3)} s")
        plt.imshow(u_plot, extent=[0, L, 0, L], origin='lower', cmap='hot')
        plt.colorbar(label='Drug Concentration')
        plt.xlabel('x (cm)')
        plt.ylabel('y (cm)')
        plt.tight_layout()
        plt.savefig(os.path.join(photo_dir, f"diffusion_{i}.png"))
        plt.close()

    fig, ax = plt.subplots(figsize=(5, 4))
    cax = ax.imshow(results[0], extent=[0, L, 0, L], origin='lower', cmap='hot')
    fig.colorbar(cax, label='Drug Concentration')
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('y (cm)')

    def update(frame):
        cax.set_array(results[frame])
        ax.set_title(f"{label} - t = {round(frame * dt * (Nt // 10), 3)} s")
        return cax,

    ani = animation.FuncAnimation(fig, update, frames=len(results), blit=False)
    ani.save(os.path.join(output_dir, f"{label}_animation.gif"), writer='pillow', fps=2)
    plt.close()

run_simulation(0.01, 0.01, 0.0, "isotropic_no_reaction")
run_simulation(0.01, 0.01, 0.1, "isotropic_with_reaction")
run_simulation(0.01, 0.001, 0.0, "anisotropic_no_reaction")
run_simulation(0.01, 0.001, 0.1, "anisotropic_with_reaction")