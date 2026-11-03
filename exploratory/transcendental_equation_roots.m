% Mathematical exploration of cot(z)+W/(z^2-1)=0 for FIXED W.
% This equation is not asserted to describe finite-well bound states.
W = 2; zmax = 10;
fun = @(z) cot(z)+W./(z.^2-1);
singularities = unique(sort([0,1,pi:pi:zmax,zmax])); roots_found = [];
for segment = 1:numel(singularities)-1
    a = singularities(segment); b = singularities(segment+1);
    margin = 1e-7*max(1,b); grid = linspace(a+margin,b-margin,2000);
    values = fun(grid);
    for k = 1:numel(grid)-1
        if isfinite(values(k)) && isfinite(values(k+1)) && values(k)*values(k+1)<0
            root = fzero(fun,[grid(k),grid(k+1)]);
            if abs(fun(root))<1e-7 && (isempty(roots_found)||all(abs(roots_found-root)>1e-6))
                roots_found(end+1) = root; %#ok<AGROW>
            end
        end
    end
end
fprintf('Roots for fixed W=%.3f in (0,%.3f):\n',W,zmax); disp(roots_found);
% A sign-change scan need not detect tangential/even-multiplicity roots.
figure; hold on;
for segment = 1:numel(singularities)-1
    grid = linspace(singularities(segment)+1e-4,singularities(segment+1)-1e-4,1000);
    plot(grid,fun(grid),'b');
end
plot(roots_found,zeros(size(roots_found)),'ro'); ylim([-10,10]);
xlabel('z'); ylabel('cot(z)+W/(z^2-1)'); title('Roots on continuous intervals'); grid on;
