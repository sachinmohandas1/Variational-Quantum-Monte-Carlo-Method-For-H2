# Variational Quantum Monte Carlo: Hydrogen Molecule Ground State Estimation

## Overview

This project uses a Variational Quantum Monte Carlo (VQMC) method to estimate the ground state energy and binding energy of a hydrogen molecule (H₂) in atomic units. It was developed as part of a computational physics course (PHYS 7411), and demonstrates the use of quantum trial wavefunctions, the Metropolis algorithm, and parameter optimization via Monte Carlo sampling.

The project uses a simplified, unitless Hamiltonian derived from atomic units, with fixed proton positions and mobile electrons. The method builds an energy function dependent on the inter-electron distance, incorporating trial wavefunctions with variational parameters to approximate the system’s ground state.

---

## Key Concepts

- **Variational principle**: Uses trial wavefunctions to estimate upper bounds on the true ground state energy.
- **Quantum Monte Carlo**: Employs the Metropolis algorithm for probabilistic sampling of quantum states.
- **Binding energy estimation**: Converts unitless energy to electronvolts to compare with experimental values.
- **Parameter tuning**:
  - `a` (scaling parameter in hydrogenic orbital): ~0.84
  - `β` (screening parameter in correlation function): optimized via energy minimization
  - `δ` (Metropolis step size): tuned to ~0.6 for ~50% acceptance rate

---

## Methodology

1. **Physical Setup**:
   - Two protons placed at fixed positions ±s/2 along the x-axis (s = 0.74 Å / a₀ ≈ 1.4 unitless).
   - Two electrons initialized with random positions in 3D space.
   - Units are normalized to the Bohr radius (a₀) and Rydberg energy.

2. **Wavefunction Construction**:
   - Trial wavefunction:
     - ψ(r₁, r₂) = φ₁(r₁) · φ₂(r₂) · f(r₁, r₂)
     - φ₁, φ₂: electron-proton orbitals (linear combinations of hydrogenic orbitals)
     - f: inter-electron correlation factor:  
       \[
       f(r_{12}) = \exp\left(\frac{r_{12}}{\alpha(1 + \beta r_{12})}\right)
       \]

3. **Metropolis Sampling**:
   - Each coordinate (x, y, z) of each electron is updated using the Metropolis acceptance criterion.
   - Sampling is repeated across a range of β values to minimize energy.

4. **Energy Evaluation**:
   - The local energy is computed analytically from the Hamiltonian and wavefunction.
   - Final energy in eV is calculated using:
     \[
     E_{\text{final}} = (\epsilon + \frac{1}{s}) \cdot 2R_y
     \]
   - Where ε is the unitless energy, s is proton separation, and R_y is the Rydberg constant in eV.

---

## Results

- **Optimized β**: ~1.1  
- **Estimated ground-state energy**: Ranged between **-30 eV to -35 eV**, depending on trial  
- **Estimated binding energy**: Centered around **4.5 eV ± 0.001–0.005 eV**, comparable to experimental H₂ value  
- **Acceptance rate**: 40–60% for δ = 0.6  
- **Stabilization strategy**: Only every 11th accepted configuration was retained to reduce autocorrelation.

A typical plot of unitless energy vs. β shows a distinct minimum, indicating the optimal variational parameter.

---

## Author

**Sachin Mohandas**  
Louisiana State University – PHYS7411: Computational Physics  
March 2022  
