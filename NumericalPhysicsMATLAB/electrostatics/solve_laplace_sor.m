function [potential,info] = solve_laplace_sor(initial,dx,dy,omega,tol,max_iter)
%SOLVE_LAPLACE_SOR Solve Laplace's equation with fixed Dirichlet boundaries.
% Rows represent y; columns represent x. tol bounds the discrete Laplacian
% residual in potential/length^2 units, not the iterate-to-iterate change.
validateattributes(initial,{'numeric'},{'2d','real','finite','nonempty'});
validateattributes(dx,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(dy,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(omega,{'numeric'},{'scalar','real','>',0,'<',2});
validateattributes(tol,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(max_iter,{'numeric'},{'scalar','integer','positive'});
[ny,nx] = size(initial);
if min(ny,nx)<3, error('Grid must have at least three rows and columns.'); end
potential = double(initial); ax = 1/dx^2; ay = 1/dy^2;
converged = false;
for iteration = 1:max_iter
    for iy = 2:ny-1
        for ix = 2:nx-1
            target = (ax*(potential(iy,ix-1)+potential(iy,ix+1)) + ...
                ay*(potential(iy-1,ix)+potential(iy+1,ix)))/(2*(ax+ay));
            potential(iy,ix) = potential(iy,ix)+omega*(target-potential(iy,ix));
        end
    end
    lap = ax*(potential(2:end-1,1:end-2)-2*potential(2:end-1,2:end-1)+ ...
        potential(2:end-1,3:end)) + ay*(potential(1:end-2,2:end-1)- ...
        2*potential(2:end-1,2:end-1)+potential(3:end,2:end-1));
    residual = max(abs(lap(:)));
    if residual <= tol, converged = true; break; end
end
info = struct('iterations',iteration,'residual',residual,'converged',converged);
if ~converged, warning('SOR did not reach residual tolerance in %d iterations.',max_iter); end
end
