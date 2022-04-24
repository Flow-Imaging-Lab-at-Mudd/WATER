function [fres, dI, dI_box, dI_gss, dI0, bias_box, bias_gss, di, di_box, ...
    di_gss, di0, mag_bias_box, mag_bias_gss, di_sd, di_sd_box, di_sd_gss] = ...
    impulse_resol(l, vr, u0, min_fres, max_fres, fres_inc, origin, props, err_level, ...
        num_ite, window_params, show_plot)
% Variation of error with feature resolution. Parameters held constant are
% 'fr', 'origin', and 'props'. At low resolutions, depending on whether the
% spacing perfectly divides the (-1, 1) region, e.g. s = 0.1 does, s = 0.3
% does not, rather discrepant behavior can be observed.
%
% Note that the proportions of error given in 'props' must include 0 as the
% first entry, so that the baseline resolution errors can be computed. When
% the error over various noise levels are averaged, supposing more than one
% nonzero noise level is given in props, the average does not include the
% baseline error with no noise.
% 
% Derek Li, June 2021

% Radius of vortex.
fr = l*vr;
% Desired feature resolutions.
fres = min_fres: fres_inc: max_fres;

% Configure windowing parameters, if applicable.
windowing = false;
if isvector(window_params) && length(window_params) == 2
    windowing = true;
    winsize = window_params(1);
    overlap = window_params(2);
end

% Generate range of spacings for evenly spaced feature resolutions.
if ~windowing
    % Global spacing of downsampled data.
    sps = zeros(1, floor((max_fres-1)/fres_inc) + 1);
    % Minimal feature resolution.
    sps(1) = fr / min_fres;
    for i = 2: size(sps, 2)
        sps(i) = fr / (fres_inc + fr/sps(i-1));
    end
else
    % Compute initial resolutions to produce the desired given feature
    % resolutions after downsampling.
    op = round(winsize .* overlap);
    op = min([op; winsize-1], [], 1);
    sps = 2*fr ./ (2*fres * (winsize-op) + winsize - 1);
end
sps_count = size(sps, 2);
% Normalize spacing.
spr = sps / fr;
% Record feature resolutions.
fres = zeros(1, sps_count);

% Vortex parameters.
density = 1000;
len_unit = 1e-3;
% Theoretical impulse values.
I0 = Hill_Impulse(density, len_unit, fr, u0);

% Containers of error.
dI = zeros(3, sps_count);
dI_box = zeros(3, sps_count);
dI_gss = zeros(3, sps_count);

dI_sd = zeros(3, sps_count);
dI_sd_box = zeros(3, sps_count);
dI_sd_gss = zeros(3, sps_count);

% Magnitudes.
di = zeros(1, sps_count);
di_box = zeros(1, sps_count);
di_gss = zeros(1, sps_count);
di_sd = zeros(1, sps_count);
di_sd_box = zeros(1, sps_count);
di_sd_gss = zeros(1, sps_count);

di0 = zeros(1, sps_count);
mag_bias_box = zeros(1, sps_count);
mag_bias_gss = zeros(1, sps_count);

% Containers for data across all runs.
% Errors here are mean absolute errors.
err = zeros(3, sps_count, num_ite);
err_box = zeros(3, sps_count, num_ite);
err_gss = zeros(3, sps_count, num_ite);

bias_box = zeros(3, sps_count, 1);
bias_gss = zeros(3, sps_count, 1);

% Error due purely to the imperfection of resolution and origin selection.
dI0 = zeros(3, sps_count);

for k = 1: sps_count
    % Construct Hill vortex with specified resolution.
    [x, y, z, u, v, w] = Hill_Vortex(spr(k), l, vr, u0, 1);
    vf = VelocityField.importCmps(x, y, z, u, v, w);
    if windowing
        vf = vf.downsample(winsize, overlap, 0);
    end
    % Focus on vortical region.
    vf.setRangePosition(fr*repmat([-1 1], 3, 1))
    fres(k) = (vf.span(1)-1)/2;
    
    [dI(:,k), dI_box(:,k), dI_gss(:,k), dI0(:,k), bias_box(:,k), bias_gss(:,k), ...
        dI_sd(:,k), dI_sd_box(:,k), dI_sd_gss(:,k), di(:,k), di_box(:,k), di_gss(:,k), ...
        di0(k), mag_bias_box(k), mag_bias_gss(k), di_sd(k), di_sd_box(k), di_sd_gss(k)] = ...
        impulse_err_run_constN(vf, props, origin, I0, num_ite, window_params, false);
end

% Compute minimum resolution needed when smoothers are applied. First row
% corresponds to biases; second row, noisy error. Box, first column;
% Gaussian, second.
min_res = -1*ones(2, 2);

% Lowest feature resolutions to achieve the desired error level.
try
    min_res(1,1) = fres(find(di0 < err_level, 1));
catch
end
try
    min_res(1,2) = fres(find(mag_bias_gss < err_level, 1));
catch
end

% try
%     min_res(2,1) = fres(find(di_box < err_level, 1));
% catch
% end
try
    min_res(2,2) = fres(find(di_gss < err_level, 1));
catch
end

disp('---------Impulse minimal resolutions---------')
% Print resolution requirements.
fprintf('Minimum resolution for %.0f%% bias: \n', err_level*100)
fprintf('Original: %d \n', min_res(1,1))
fprintf('Gaussian: %d \n', min_res(1,2))

fprintf('Minimum resolution for %.0f%% noise-propagated error: \n', err_level*100)
% fprintf('Box: %d \n', min_res(2,1))
fprintf('Gaussian: %d \n', min_res(2,2))

%%%%%%%%%%%%%%%% Dimensional Plots %%%%%%%%%%%%%%%%%
% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [];
dim_str = {'x', 'y', 'z'};

if ~exist('show_plot', 'var') || ~show_plot
    return
end

for dim = dims
    % Smoother bias plot.
    figure;
    scatter(fres, dI0(dim,:), 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    scatter(fres, bias_box(dim,:), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    scatter(fres, bias_gss(dim,:), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
    
    legend({'imperfect resolution', ...
        'box-filtered', ...
        'Gaussian-filtered'}, ...
        'Interpreter', 'latex')
    xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
    ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
    title(strcat('$', dim_str{dim}, '$ Smoother Bias at $r = $', {' '}, string(fr)))
    
    % Mean error plot.
    figure;
    errorbar(fres, dI(dim,:), dI_sd(dim,:), 'ko', 'MarkerFaceColor','black', 'LineWidth', 1)
    hold on
    errorbar(fres, dI_box(dim,:), dI_box_sd(dim,:), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
    hold on
    errorbar(fres, dI_gss(dim,:), dI_gss_sd(dim,:), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)
    
    legend({'unfiltered', ...
        'box-filtered', ...
        'Gaussian-filtered'}, ...
        'Interpreter', 'latex')
    xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
    ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
    title(strcat('$', string(dim_str{dim}), '$ Mean Error at $\delta u = $', ...
        string(props(2)*100), '\% at $r = $', {' '}, string(fr)))
end


%%%%%%%%%%%%%%%%%%% Magnitude Plots %%%%%%%%%%%%%%%%%%%%

% Smoother bias plot.
figure;
scatter(fres, di0, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
hold on
scatter(fres, mag_bias_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
scatter(fres, mag_bias_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
hold on

legend({'box filtered', ...
    'Gaussian-filtered'}, ...  
    'Interpreter', 'latex')
xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
ylabel('$\left|\frac{\delta I}{I}\right|$')
title(strcat('Magnitude of Smoother Bias at $r = $', {' '}, string(fr)))

% Mean error plot.
figure;
errorbar(fres, di, di_sd, 'ko', 'MarkerFaceColor','black', 'LineWidth', 1)
hold on
errorbar(fres, di_box, di_box_sd, 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
errorbar(fres, di_gss, di_gss_sd, 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)
hold on

legend({'unfiltered', ...
    'box-filtered', ...
    'Gaussian-filtered'}, ...  
    'Interpreter', 'latex')
xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
ylabel('$\left|\frac{\delta I}{I}\right|$')
title(strcat('Mean Error Magnitude at $\delta u = $', ...
        string(props(2)*100), '\% at $r = $', {' '}, string(fr)))