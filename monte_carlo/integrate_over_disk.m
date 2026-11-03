% Integral of sin(log(x+y+1)) over x^2+y^2 < 1/4.
rng(42); N = 1e5; bounds = [-0.5,0.5; -0.5,0.5];
points = bounds(:,1)' + rand(N,2).*diff(bounds,1,2)';
inside = sum(points.^2,2) < 0.25;
values = zeros(N,1);
% Evaluate only inside the disk, where x+y+1 is strictly positive.
values(inside) = sin(log(sum(points(inside,:),2)+1));
area = prod(diff(bounds,1,2));
integral_estimate = area*mean(values);
integral_se = area*std(values,0)/sqrt(N);
fprintf('Integral: %.8f +/- %.8f (one sampling SE)\n',integral_estimate,integral_se);
