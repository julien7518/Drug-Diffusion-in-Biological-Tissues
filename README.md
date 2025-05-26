


# Drug Diffusion in Biological Tissues

This project models the diffusion of a drug through biological tissue using both analytical and numerical methods. It is implemented in Python and simulates 2D diffusion under different conditions (isotropic, anisotropic, with and without reaction).

## Contents

- `report.pdf`: Final academic report describing the theory, methods, and results.
- `drug_diffusion_simulation.py`: Main simulation script.
- `generate_analytical_solution.py`: Script to generate the analytical solution figure.
- `resources/`: Contains output images and animations from the simulations.

## Features

- Analytical solution of the 2D diffusion equation.
- Finite difference numerical simulation (explicit Euler).
- Comparison between isotropic and anisotropic diffusion.
- Option to include a reaction term.
- Visual outputs: PNG heatmaps and animated GIFs.
- Convergence test and RMSE evaluation.

## Usage

Run the simulation with:

```bash
python drug_diffusion_simulation.py
```

Generate the analytical figure with:

```bash
python generate_analytical_solution.py
```