% Basic MATLAB numerical checks. Run from the repository root.
project_root = fileparts(mfilename('fullpath'));
addpath(fullfile(project_root,'electrostatics'),fullfile(project_root,'quantum_wells'));
% Square-well roots: positive, ordered, below barrier, satisfy a parity equation.
m = 9.10938356e-31; q = 1.602176634e-19; hbar = 1.0545718e-34;
for depth_ev = [0.001,0.1,1]
    V0 = depth_ev*q; L = 1e-9;
    levels = finite_well_levels(m,V0,L,hbar,q);
    assert(~isempty(levels) && all(levels>0) && all(levels<V0/q));
    assert(all(diff(levels)>0));
    z = (L/2)*sqrt(2*m*levels*q)/hbar;
    rho = (L/2)*sqrt(2*m*V0)/hbar;
    decay = sqrt(rho^2-z.^2);
    residual = min(abs(z.*tan(z)-decay),abs(-z./tan(z)-decay));
    assert(max(residual)<1e-6);
    assert(numel(levels)==ceil(2*rho/pi));
end
% Unequal grid spacing, fixed boundaries and known harmonic solution.
x = linspace(0,1,17); y = linspace(0,1,13); [X,Y] = meshgrid(x,y);
exact = X.*Y; initial = zeros(size(exact));
initial([1,end],:) = exact([1,end],:); initial(:,[1,end]) = exact(:,[1,end]);
[answer,info] = solve_laplace_sor(initial,x(2)-x(1),y(2)-y(1),1.5,1e-8,20000);
assert(info.converged && info.residual<=1e-8);
assert(max(abs(answer(:)-exact(:)))<1e-7);
assert(isequal(answer([1,end],:),initial([1,end],:)));
assert(isequal(answer(:,[1,end]),initial(:,[1,end])));
% Sine/sinh solution: discretisation error decreases on a refined grid.
errors = zeros(1,2);
for trial = 1:2
    n = 10*2^(trial-1)+1; x = linspace(0,1,n); [X,Y] = meshgrid(x,x);
    exact = sin(X).*sinh(Y); initial = zeros(n);
    initial([1,end],:) = exact([1,end],:); initial(:,[1,end]) = exact(:,[1,end]);
    answer = solve_laplace_sor(initial,x(2)-x(1),x(2)-x(1),1.6,1e-8,20000);
    errors(trial) = max(abs(answer(:)-exact(:)));
end
assert(errors(2)<errors(1)/3);
fprintf('Finite-well and Laplace validation checks passed.\n');
