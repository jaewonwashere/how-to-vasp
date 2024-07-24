# VASP

A beginner's guide to adsorption energy calculations using VASP.

## Input Files

There are four main input files for you to run a VASP calculation: INCAR, POSCAR, POTCAR, and KPOINTS.

Let's get right into it.

### KPOINTS

When we use DFT to calculate the electron density, we need to integrate over space. For periodic systems (which we can construct with supercells), the solutions of the Schrodinger equation satisfy [Bloch's theorem](https://en.wikipedia.org/wiki/Bloch%27s_theorem).

* Bloch function (the form that a wave function in a periodic potential takes):

  ![image](https://github.com/user-attachments/assets/6e0966cb-0ffe-4b77-bb4f-5ad730a563d5)

This property lets us solve the Schrodinger equation for each value of the wave vector, **k**.

To do this, we need to move from the 3d-space that we live in, called real space, to the space where the wave vectors live. This space is called the *k-space* (or the reciprocal space). The term *k-space* makes sense: it is where the wave vectors, **k** live. Why is it called the reciprocal space? I won't go into the details here, but here's the takeaway: if you send a vector, **a**, from real space to reciprocal space, its magnitude changes from |**a**| to 2π/|**a**|. Hence, the name 'reciprocal'. (If you want to learn more, Chapter 3.1 in Sholl and Steckel's book has a nice explanation).

**So, why do we want to use the k-space? The computer likes it!**

If we were to integrate in real space, we would have to do a continuous integral, because the space we are integrating over is continuous. In reciprocal space, the integral is only calculated over possible values of **k**. This makes the calculation a lot less intense.

For the computer to run these calculations, we need to tell it how to do it. This is done by specifying a *k-point mesh*. The most common way to do this is using a Monkhorst-Pack grid (essentially a uniformly spaced grid in the [Brillouin Zone](https://en.wikipedia.org/wiki/Brillouin_zone)). Using this method, if we tell the computer the number of *k-points* to sample in each direction (_N<sub>1</sub>_ × _N<sub>2</sub>_ × _N<sub>3</sub>_). These intergers are called 'subdivisions'. VASP (or other DFT package) will take care of the rest.

As you may have guessed, using a denser *k-point grid* will lead to more accurate, but more computationally intensive results. A tradeoff! Also, remember that your cell in the reciprocal space scales as the inverse of your unit cell in real space. This means that having a larger unit cell in real space will require less *k-points*, and vice versa.

Thus, we need to do a **convergence test** to see which value of *n* would give us a well-converged calculation for our given structure.

The VASP input file for _k-points_ is KPOINTS. Let's look at what this looks like.

```
Regular k-point mesh
0              ! 0 -> determine number of k points automatically
Monkhorst-Pack ! Also works with M or m (other option: Gamma for Gamma centered mesh (G, g))
4  4  4        ! subdivisions N_1, N_2 and N_3 along the reciprocal lattice vectors
0  0  0        ! optional shift of the mesh (s_1, s_2, s_3)
```

