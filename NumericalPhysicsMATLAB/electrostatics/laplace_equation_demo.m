% Rectangular-grid demonstrations with consistent meshgrid orientation.
nx = 41; ny = 31; lx = 1; ly = 1;
x = linspace(0,lx,nx); y = linspace(0,ly,ny);
dx = x(2)-x(1); dy = y(2)-y(1); [X,Y] = meshgrid(x,y);
exact_fields = {sin(X).*sinh(Y), X.*Y};
for example = 1:2
    exact = exact_fields{example}; initial = zeros(ny,nx);
    initial([1,end],:) = exact([1,end],:);
    initial(:,[1,end]) = exact(:,[1,end]);
    [potential,info] = solve_laplace_sor(initial,dx,dy,1.7,1e-8,20000);
    fprintf('Example %d: %d iterations; residual %.3g; max error %.3g\n', ...
        example,info.iterations,info.residual,max(abs(potential(:)-exact(:))));
    figure;
    subplot(1,2,1); contourf(X,Y,potential,30); colorbar;
    title('Numerical potential'); xlabel('x'); ylabel('y');
    subplot(1,2,2); contourf(X,Y,exact,30); colorbar;
    title('Analytical potential'); xlabel('x'); ylabel('y');
end
