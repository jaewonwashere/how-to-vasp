# VASP: A Beginner's Guide to Adsorption Energy Calculations

This guide introduces the basics of performing adsorption energy calculations using VASP (Vienna Ab initio Simulation Package).

## Table of Contents
1. [Input Files](#input-files)
2. [KPOINTS](#kpoints)
   - [Understanding K-space](#understanding-k-space)
   - [K-point Mesh](#k-point-mesh)
   - [Convergence Testing](#convergence-testing)

## Input Files

To run a VASP calculation, you need four main input files:

1. INCAR
2. POSCAR
3. POTCAR
4. KPOINTS

Let's explore these files, starting with KPOINTS.

## KPOINTS

### Understanding K-space

In Density Functional Theory (DFT) calculations for periodic systems, we utilize k-space (reciprocal space) instead of real space. This approach is based on Bloch's theorem and offers computational advantages.

#### Bloch's Theorem
The wave function in a periodic potential takes the form of a Bloch function:

![Bloch Function](https://github.com/user-attachments/assets/6e0966cb-0ffe-4b77-bb4f-5ad730a563d5)

Key points about k-space:
- It's where wave vectors (**k**) reside
- Also known as reciprocal space
- Vectors in k-space are inversely proportional to real space vectors
  - |**a**| in real space → 2π/|**a**| in k-space

> 💡 **Tip:** For a deeper understanding of reciprocal space, refer to Chapter 3.1 in Sholl and Steckel's book.

### K-point Mesh

To perform calculations in k-space, we define a k-point mesh. The most common method is the Monkhorst-Pack grid, which creates a uniformly spaced grid in the [Brillouin Zone](https://en.wikipedia.org/wiki/Brillouin_zone).

#### Specifying the Mesh
- Define the number of k-points in each direction: _n<sub>1</sub>_ × _n<sub>2</sub>_ × _n<sub>3</sub>_
- VASP handles the rest of the setup

#### Trade-offs
- Denser k-point grid → More accurate but computationally intensive
- Larger real-space unit cell → Fewer k-points needed (and vice versa)

### Convergence Testing

To determine the optimal k-point mesh density, perform a convergence test:
1. Run calculations with increasing numbers of k-points
2. Analyze the results to find the point where increasing k-points no longer significantly improves accuracy
3. Choose the k-point density that balances accuracy and computational cost

> ⚠️ **Important:** Always perform a convergence test for your specific structure to ensure reliable results.
