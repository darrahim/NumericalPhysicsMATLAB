# Independent numerical reference checks

These results were obtained using Python/NumPy/SciPy implementations of the revised equations. They are not executions of the MATLAB files. MATLAB syntax and toolbox-based fitting remain to be checked in MATLAB.

- Well depth 0.001 eV: 1 bound state(s), energies [0.000993495] eV; maximum dimensionless residual 2.63e-15.
- Well depth 0.1 eV: 1 bound state(s), energies [0.063707032] eV; maximum dimensionless residual 1.11e-16.
- Well depth 1 eV: 2 bound state(s), energies [0.190760543, 0.702769568] eV; maximum dimensionless residual 1.11e-15.
- Rectangular Laplace xy case: maximum error 2.6e-10; residual 9.15e-09.
- Sin(x)sinh(y) grid refinement: errors 3.46e-05 -> 8.77e-06; reduction factor 3.947.
- Full-basis wavefunction evolution: maximum discrete norm error 7.77e-16.
- Triangular-aperture transform: maximum normalised quadrature error 8.33e-08; analytical first zero 4*pi/w.

Run `run_validation.m` in MATLAB to check the actual finite-well and Laplace implementations. The Monte Carlo scripts report their own sampling errors; no external spectroscopy dataset was supplied, so the spectral fitter has not been validated on experimental data.
