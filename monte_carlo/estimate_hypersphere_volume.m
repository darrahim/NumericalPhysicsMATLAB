% Estimate a D-dimensional unit-ball volume and its sampling standard error.
rng(42); D = 10; N = 1e6; lower = -1; upper = 1;
points = lower + (upper-lower)*rand(N,D);
inside = sum(points.^2,2) <= 1;
box_volume = (upper-lower)^D;
volume_estimate = box_volume*mean(inside);
volume_se = box_volume*std(double(inside),0)/sqrt(N);
volume_exact = pi^(D/2)/gamma(D/2+1);
fprintf('Volume: %.6f +/- %.6f (one sampling SE); exact: %.6f\n', ...
    volume_estimate,volume_se,volume_exact);
% A zero acceptance count gives an uninformative zero plug-in SE.
if ~any(inside) || all(inside)
    warning('Degenerate acceptance count: increase N before interpreting uncertainty.');
end
