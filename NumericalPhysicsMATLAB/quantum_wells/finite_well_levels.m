function E_levels = finite_well_levels(m,V0,L,hbar,q)
%FINITE_WELL_LEVELS Bound energies in eV, measured from the well bottom.
% Symmetric well of full width L: V=0 inside and V=V0 outside.
% Equal effective mass inside/outside. Inputs: kg, J, m, J*s, J/eV.
% With z=kL/2 and rho=(L/2)*sqrt(2*m*V0)/hbar, the branches are
% even: z*tan(z)=sqrt(rho^2-z^2); odd: -z*cot(z)=sqrt(rho^2-z^2).
% Each continuous positive branch is bracketed separately; no pole is crossed.
validateattributes(m,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(V0,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(L,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(hbar,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(q,{'numeric'},{'scalar','real','finite','positive'});
rho = (L/2)*sqrt(2*m*V0)/hbar; roots = [];
for branch = 0:ceil(2*rho/pi)-1
    left = branch*pi/2; right = min((branch+1)*pi/2,rho);
    if right <= left, continue; end
    margin = min((right-left)/4,32*eps(max(1,right)));
    lo = left+margin; hi = right-margin;
    if mod(branch,2)==0
        equation = @(z) z.*tan(z)-sqrt(max(0,rho^2-z.^2));
    else
        equation = @(z) -z./tan(z)-sqrt(max(0,rho^2-z.^2));
    end
    if equation(lo)<0 && equation(hi)>0
        z = fzero(equation,[lo,hi],optimset('TolX',1e-12));
        if abs(equation(z))>1e-7*max(1,rho)
            error('Root residual too large; inspect numerical conditioning.');
        end
        roots(end+1) = z; %#ok<AGROW>
    else
        warning('Unresolved branch near a threshold; check extremely shallow states.');
    end
end
E_levels = sort(V0*(roots/rho).^2/q);
end
