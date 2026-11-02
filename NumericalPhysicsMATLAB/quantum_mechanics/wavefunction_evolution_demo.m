% Dirichlet walls at x=0,L; all Hamiltonian degrees of freedom are interior.
L = 1; m = 1; hbar = 1; N = 500;
dx = L/(N+1); x = (1:N)'*dx;
potential = sin(4*pi*x/L).^4;
kinetic = -(hbar^2/(2*m*dx^2))*spdiags([ones(N,1),-2*ones(N,1),ones(N,1)],-1:1,N,N);
H = kinetic+spdiags(potential,0,N,N);
[eigenfunctions,eigenvalues_matrix] = eig(full(H));
[energies,order] = sort(diag(eigenvalues_matrix));
eigenfunctions = eigenfunctions(:,order)/sqrt(dx);
psi_initial = sqrt(2/L)*sin(pi*x/L);
coefficients = eigenfunctions' * psi_initial * dx;
% Full eigenbasis preserves the discrete norm (no unquantified truncation).
times = linspace(0,5,6); states = eigenfunctions*(coefficients.*exp(-1i*energies*times/hbar));
fprintf('Max discrete norm error: %.3g\n',max(abs(dx*sum(abs(states).^2,1)-1)));
figure; plot([0;x;L],[zeros(1,5);eigenfunctions(:,1:5);zeros(1,5)]);
xlabel('x'); ylabel('Eigenfunction'); title('Lowest five energy eigenfunctions');
figure; plot([0;x;L],[zeros(1,numel(times));abs(states).^2;zeros(1,numel(times))]);
xlabel('x'); ylabel('|psi(x,t)|^2'); title('Probability density evolution');
legend(arrayfun(@(t) sprintf('t=%.2f',t),times,'UniformOutput',false));
% Energies and wavefunctions are discretised approximations: check grid refinement.
