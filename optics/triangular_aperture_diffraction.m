% Single triangular aperture with FULL support width w.
% A(y)=max(1-|y|/a,0), a=w/2.
% Integral A(y)exp(-iky)dy = a*[sin(ka/2)/(ka/2)]^2.
w = 20e-6; a = w/2; y = linspace(-a,a,4001);
aperture = max(1-abs(y)/a,0);
k = linspace(-12*pi/w,12*pi/w,601);
u = k*a/2; sinc_unnormalised = ones(size(u));
nonzero = abs(u)>eps; sinc_unnormalised(nonzero)=sin(u(nonzero))./u(nonzero);
amplitude_analytic = a*sinc_unnormalised.^2;
amplitude_numeric = zeros(size(k));
for index = 1:numel(k)
    amplitude_numeric(index) = real(trapz(y,aperture.*exp(-1i*k(index)*y)));
end
fprintf('First positive zero: k=4*pi/w = %.6g rad/m\n',4*pi/w);
fprintf('Max normalised transform error: %.3g\n',max(abs(amplitude_numeric-amplitude_analytic))/a);
figure; subplot(2,1,1); plot(y,aperture); xlabel('y (m)'); ylabel('Aperture amplitude');
subplot(2,1,2); plot(k,amplitude_analytic/a,k,amplitude_numeric/a,'--');
xlabel('Spatial angular frequency k (rad/m)'); ylabel('Normalised field amplitude');
legend('Analytical','Numerical'); title('Single triangular aperture Fourier transform');
