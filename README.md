# Numerical Physics in MATLAB

**Darius Rahimifard**

Numerical physics projects exploring Monte Carlo integration, electrostatics, quantum confinement, wavefunction evolution and diffraction, alongside a tool for fitting experimental photoluminescence spectra.

## Projects

| Folder | Entry point | Method and output |
| --- | --- | --- |
| `monte_carlo` | `estimate_hypersphere_volume.m` | Uniform sampling of a 10D unit ball; estimate, sampling standard error and exact volume comparison. |
| `monte_carlo` | `estimate_toroidal_mass_properties.m` | Mass, centre of mass and origin-axis moments for a restricted toroidal region, including sampling uncertainties. |
| `monte_carlo` | `integrate_over_disk.m` | Two-dimensional integral with a circular domain and sampling error. |
| `electrostatics` | `laplace_equation_demo.m` | SOR on rectangular grids; discrete residual, analytical comparison and contour plots. |
| `quantum_wells` | `predict_quantum_well_pl.m` | Approximate GaAs/AlGaAs transition energies from finite-well confinement and an empirical exciton correction. |
| `quantum_wells` | `fit_photoluminescence_spectra.m` | Joint one/two-Gaussian fitting with constant background, residuals and approximate confidence intervals. Requires experimental data. |
| `quantum_mechanics` | `wavefunction_evolution_demo.m` | Finite-difference Hamiltonian with Dirichlet walls and evolution in its full eigenbasis. |
| `optics` | `triangular_aperture_diffraction.m` | Analytical triangular-aperture Fourier transform compared with numerical quadrature. |
| `exploratory` | `transcendental_equation_roots.m` | Numerical roots of a fixed mathematical equation on intervals separated by its singularities. |

## Quick start

Use MATLAB R2019b or newer. Open the repository folder in MATLAB and run:

```matlab
run_validation
```

For a demonstration, change the Current Folder to the relevant project folder and run its entry-point script. For example:

```matlab
cd electrostatics
laplace_equation_demo
```

Keep each project folder's files together. `solve_laplace_sor.m` and `finite_well_levels.m` are reusable functions; their companion demonstrations show how to call them. The Monte Carlo scripts set `rng(42)` for repeatable samples. The hypersphere example defaults to one million samples: reduce `N` for a quick run.

The spectral fitter needs Curve Fitting Toolbox and Signal Processing Toolbox. It asks for a folder of `.asc` files, each containing two numeric columns: pixel coordinate and counts. Those experimental data are not bundled. All other projects use base MATLAB.

## Numerical conventions and limitations

- Finite-well energies are measured from the well bottom, with a barrier V0 and full well width L. The model uses equal mass inside and outside the well. Inputs are SI; output energies are eV. Roots are bracketed on continuous parity branches and checked by residual. Extremely shallow states close to machine precision can remain unresolved and produce a warning.
- The quantum-well PL model uses approximate material parameters and an empirical exciton interpolation; agreement with experiment is not established here.
- Laplace matrices use rows for y and columns for x. Tolerance controls the discrete Laplacian residual; continuum accuracy also depends on grid spacing.
- Monte Carlo errors are estimated one-standard-error sampling uncertainties. Centre-of-mass errors use the delta method, including numerator/denominator covariance. The toroidal region is restricted by the script's stated box; inertia is about the origin. Degenerate acceptance counts require more samples.
- Wavefunction evolution conserves the discrete norm using the full basis; continuum energies require grid-convergence checks.
- The spectral fitter selects at most two prominent peaks. More complex spectra, nonconstant backgrounds and counting-noise weighting may need a different model. Confidence intervals rely on least-squares assumptions; pixel-to-energy conversion requires experimental calibration.
- The optics example is a single triangular aperture. The first nonzero transform zero is 4*pi/w when w is the full support width.
- The exploratory transcendental equation is not presented as a derived quantum-well model. Its sign-change search does not guarantee detection of tangential roots.


Original coursework was reorganised and revised with AI assistance. 
