% Uniform-density toroidal region restricted to x in [1,4], y in [-3,sqrt(15)].
% Moments are about the coordinate origin, not the centre of mass.
rng(42); N = 1e5; density = 1;
bounds = [1,4; -3,sqrt(15); -1,1];
points = bounds(:,1)' + rand(N,3).*diff(bounds,1,2)';
box_volume = prod(diff(bounds,1,2));
inside = points(:,3).^2 + (hypot(points(:,1),points(:,2))-3).^2 <= 1;
weight = density*double(inside);
mass = box_volume*mean(weight);
mass_se = box_volume*std(weight,0)/sqrt(N);
if mass == 0, error('No accepted points; increase N or check bounds.'); end
centre = sum(points.*weight,1)/sum(weight);
% Delta-method influence values include numerator/denominator covariance.
influence = weight.*(points-centre)/mean(weight);
centre_se = std(influence,0,1)/sqrt(N);
squared_distances = [points(:,2).^2+points(:,3).^2, ...
    points(:,1).^2+points(:,3).^2, points(:,1).^2+points(:,2).^2];
integrands = squared_distances.*weight;
inertia_origin = box_volume*mean(integrands,1);
inertia_se = box_volume*std(integrands,0,1)/sqrt(N);
fprintf('Mass: %.6f +/- %.6f\n',mass,mass_se);
fprintf('Centre [x y z]: %s\n',mat2str(centre,6));
fprintf('Centre SE:     %s\n',mat2str(centre_se,6));
fprintf('Origin inertia [Ixx Iyy Izz]: %s\n',mat2str(inertia_origin,6));
fprintf('Inertia SE:                  %s\n',mat2str(inertia_se,6));
