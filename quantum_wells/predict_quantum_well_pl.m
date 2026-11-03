% Approximate constant-mass confinement model; exciton term is empirical.

%% Physical constants
hbar = 1.0545718e-34;      % J·s
m0   = 9.10938356e-31;     % kg
q    = 1.602176634e-19;    % C (J/eV)

%% Material parameters (83 K)
Eg_GaAs = 1.505;           % eV

x = 0.33;
Eg_AlGaAs = Eg_GaAs + 1.247*x;
DeltaEg = Eg_AlGaAs - Eg_GaAs;

Ve = 0.65 * DeltaEg * q;   % J
Vh = 0.35 * DeltaEg * q;   % J

me = 0.067 * m0;
mhh = 0.45 * m0;

%% Well widths (nm)
L_list_nm = [10, 5, 2.5, 1.2];

%% Exciton parameters (meV)
Ex_bulk = 4.0;   % meV
Ex_2D   = Ex_bulk * 4;  % meV
L0      = 5.0;   % nm

fprintf('GaAs/Al0.33Ga0.67As PL energies at 83 K\n\n');

for L_nm = L_list_nm
    L = L_nm * 1e-9;

    Ee = finite_well_levels(me, Ve, L, hbar, q);
    Eh = finite_well_levels(mhh, Vh, L, hbar, q);

    % Exciton binding energy (eV)
    Ex_meV = Ex_bulk + (Ex_2D - Ex_bulk) * (1 - exp(-L0/L_nm));
    Ex = Ex_meV * 1e-3;

    fprintf('Well width = %.2f nm\n', L_nm);
    fprintf('  Exciton binding = %.2f meV\n', Ex_meV);

    for n = 1:min(length(Ee), length(Eh))
        EPL = Eg_GaAs + Ee(n) + Eh(n) - Ex;
        fprintf('  e%d-hh%d : %.4f eV\n', n, n, EPL);
    end
    fprintf('\n');
end
