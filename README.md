# Molecular Dynamics Simulation: Extended System (PBC & MIC)

This project is from independent research and development conducted during my **Master of Science program in Computational Physics**. The entire methodology, analysis, and implementation of all numerical solvers were executed solely by me, and the full work is documented in the accompanying report and code.

---

## Table of Contents

* [About The Project](#about-the-project)
* [Core Objectives](#core-objectives)
* [Languages and Libraries](#languages-and-libraries)
* [Methods Implemented](#methods-implemented)
* [Key Findings](#key-findings)
* [Getting Started](#getting-started)
* [Full Project Report](#full-project-report)
* [Contact](#contact)

---

## About The Project

This work extends the basic two-particle dynamics simulation into a full **Molecular Dynamics (MD) simulation** of an extended system containing multiple atoms. It utilizes the **Lennard-Jones (LJ) potential** for inter-atomic interaction and the **Velocity Verlet integrator** for time evolution.

The crucial addition is the implementation of **Periodic Boundary Conditions (PBC)** combined with the **Minimum Image Convention (MIC)**, allowing the simulation of bulk material properties by avoiding boundary effects and simulating an effectively infinite system.

## Core Objectives

1.  **Implement Periodic Boundary Conditions (PBC):** Develop an algorithm to handle particle movement across system boundaries, allowing atoms leaving one side to re-enter from the opposite.
2.  **Apply Minimum Image Convention (MIC):** Optimize the force calculation by ensuring each particle interacts only with the closest image of any other particle, including images in neighboring simulation boxes.
3.  **Establish Initial Conditions:** Implement methods for setting initial positions (e.g., a cubic lattice) and initial velocities (e.g., Maxwell-Boltzmann distribution), including **drift removal** and velocity **rescaling**.
4.  **Verify Conservation Laws:** Demonstrate that the total energy (Hamiltonian) remains conserved throughout the simulation, validating the combined PBC/Verlet approach for large systems.

---

## Languages and Libraries

| Category | Tools & Libraries | Competency Demonstrated |
| :--- | :--- | :--- |
| **Language** | Python | Efficient development of the core MD simulation loop and data handling. |
| **Numerical** | NumPy | Advanced array manipulation for handling positions, velocities, and calculating interaction matrices efficiently. |
| **Visualization** | Matplotlib | Generating high-quality plots for particle trajectories, and kinetic/potential/total energy conservation. |

---

## Methods Implemented

The computational core of the project implements the following techniques:

| Method | Role in Project | Key Implementation Detail |
| :--- | :--- | :--- |
| **Periodic Boundary Conditions (PBC)** | Boundary handling. | Enables simulation of bulk properties by making the simulation box infinitely repeating. |
| **Minimum Image Convention (MIC)** | Interaction calculation. | Used within PBC to compute the shortest distance between atoms, preventing unphysical interactions across large distances. |
| **Velocity Verlet Algorithm** | Time integration. | A stable, high-order algorithm essential for accurately integrating the equations of motion over thousands of time steps. |
| **Kinetic Energy Rescaling** | Thermal control. | Used to maintain the system at a desired average temperature (constant $\bar{K}$) by occasionally rescaling particle velocities. |

## Key Findings

* **Effective Bulk Modeling:** The combined implementation of PBC and MIC successfully allowed the simulation of continuous atomic interaction, essential for modeling material properties.
* **Energy Stability:** The total energy of the extended system showed only minor, bounded fluctuations over long time periods, confirming the robustness and stability of the Velocity Verlet integrator within the PBC environment.
* **Velocity Distribution:** The system naturally evolved from the initial, ordered state to a disordered state, and the particle velocities showed characteristics of a thermalized system.
* **Successful Trajectory Tracking:** Position and velocity plots for individual atoms confirmed physically meaningful movement within the confined, periodic space.

---

## Getting Started

### Execution

To run the simulation and generate the results and visualizations, execute the core solver script:

```bash
python Molecular_Dynamics.py
```
---

## Full Project Report

For a complete breakdown of the theoretical derivations and full numerical results, please see the final project report:

[**Full MD Extended System Project Report (PDF)**](Molecular_Dynamics.pdf)

---

## Contact

I'm happy to hear your feedback or answer any questions about this project!

**Author** Rama Khalil

**Email**  rama.khalil.990@gmail.com
