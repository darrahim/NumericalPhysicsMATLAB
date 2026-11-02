% Jointly fit one or two Gaussian peaks plus a constant background.
% Requires Curve Fitting Toolbox and Signal Processing Toolbox.
folder = uigetdir([], 'Select folder containing two-column .asc spectra');
if isequal(folder,0), fprintf('Selection cancelled.\n'); return; end
files = dir(fullfile(folder,'*.asc'));
if isempty(files), warning('No .asc files found.'); return; end
for file_index = 1:numel(files)
    data = readmatrix(fullfile(folder,files(file_index).name));
    if size(data,2)<2, warning('Skipping %s: two columns required.',files(file_index).name); continue; end
    data = data(all(isfinite(data(:,1:2)),2),1:2);
    [pixel,order] = sort(data(:,1)); counts = data(order,2);
    if numel(pixel)<10 || any(diff(pixel)<=0)
        warning('Skipping %s: need >=10 distinct pixel coordinates.',files(file_index).name); continue;
    end
    background = median(counts); smoothed = movmean(counts,7);
    prominence = 0.05*(max(smoothed)-min(smoothed));
    if prominence<=0, warning('Skipping flat spectrum.'); continue; end
    [heights,locations] = findpeaks(smoothed,pixel,'MinPeakProminence',prominence);
    if isempty(heights), warning('No peaks detected in %s.',files(file_index).name); continue; end
    [~,ranked] = sort(heights,'descend'); chosen = ranked(1:min(2,numel(ranked)));
    [centres,ordering] = sort(locations(chosen)); amplitudes = max(heights(chosen(ordering))-background,eps);
    spacing = median(diff(pixel)); span = pixel(end)-pixel(1); sigma0 = max(span/50,spacing);
    if numel(chosen)==1
        model = fittype('b+A1*exp(-0.5*((x-c1)/s1)^2)', ...
            'independent','x','coefficients',{'b','A1','c1','s1'});
        start = [background,amplitudes(1),centres(1),sigma0];
        lower = [-Inf,0,pixel(1),spacing/10]; upper = [Inf,Inf,pixel(end),span];
    else
        midpoint = mean(centres);
        model = fittype('b+A1*exp(-0.5*((x-c1)/s1)^2)+A2*exp(-0.5*((x-c2)/s2)^2)', ...
            'independent','x','coefficients',{'b','A1','c1','s1','A2','c2','s2'});
        start = [background,amplitudes(1),centres(1),sigma0,amplitudes(2),centres(2),sigma0];
        lower = [-Inf,0,pixel(1),spacing/10,0,midpoint,spacing/10];
        upper = [Inf,Inf,midpoint,span,Inf,pixel(end),span];
    end
    try
        [fitted,gof] = fit(pixel,counts,model,'StartPoint',start,'Lower',lower,'Upper',upper);
        fprintf('\n%s: R-squared %.5f; RMSE %.5f\n',files(file_index).name,gof.rsquare,gof.rmse);
        fprintf('Peak 1: centre %.4f pixels; FWHM %.4f pixels\n',fitted.c1,2*sqrt(2*log(2))*fitted.s1);
        if numel(chosen)==2
            fprintf('Peak 2: centre %.4f pixels; FWHM %.4f pixels\n',fitted.c2,2*sqrt(2*log(2))*fitted.s2);
        end
        disp('Approximate 95% coefficient confidence intervals (coefficient order follows):');
        disp(coeffnames(fitted)); disp(confint(fitted,0.95));
        prediction = feval(fitted,pixel);
        figure; subplot(2,1,1); plot(pixel,counts,'k.',pixel,prediction,'r-');
        title(files(file_index).name,'Interpreter','none'); xlabel('Pixel'); ylabel('Counts'); legend('Data','Joint fit');
        subplot(2,1,2); plot(pixel,counts-prediction); xlabel('Pixel'); ylabel('Residual');
    catch fit_error
        warning('Fit failed for %s: %s',files(file_index).name,fit_error.message);
    end
end
% Gaussian/unweighted least-squares model: confidence intervals are approximate.
% Convert to wavelength/energy only with the experimental calibration.
