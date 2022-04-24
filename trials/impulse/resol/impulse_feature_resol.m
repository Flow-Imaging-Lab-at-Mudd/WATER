% Apply impulse_resol.m at different vortical radii to generate graphs of error
% versus feature resolution, incorporating variation of radii.
%
% Derek Li, March 2022

% Constant parameters.
vr = 1;
u0 = 1;
% Windowing parameters.
window_params = [];

% Proportions of velocity noise introduced.
props = [0 1.5];
% Number of iterations per noise trial.
num_ite = 5;

% Desired level of error.
err_level = 0.1;

 % defined as 1.
min_fres = 1;
max_fres = 10;
fres_inc = 1;
% Number of uniformly spaced resolutions used.
sps_count = floor((max_fres-1)/fres_inc) + 1;

% Reference point for impulse calculation.
origin = [0 0 0]';

% Radius of the vortex to be varied in this file.
radii = 0.05: 0.1: 1;
radii_count = size(radii, 2);

% Mean errors and baseline biases.
dI = zeros(3, sps_count, radii_count);
dI_box = zeros(3, sps_count, radii_count);
dI_gss = zeros(3, sps_count, radii_count);

% Imperfection of resolution.
dI0 = zeros(3, sps_count, radii_count);
di0 = zeros(sps_count, radii_count);

bias_box = zeros(3, sps_count, radii_count);
bias_gss = zeros(3, sps_count, radii_count);
mag_bias_box = zeros(sps_count, radii_count);
mag_bias_gss = zeros(sps_count, radii_count);

di = zeros(sps_count, radii_count);
di_box = zeros(sps_count, radii_count);
di_gss = zeros(sps_count, radii_count);

% Collect data for different radii.
for i = 1: radii_count
    l = radii(i);
    % Vary global resolution for the given radius.
    [fres, dI(:,:,i), dI_box(:,:,i), dI_gss(:,:,i), dI0(:,:,i), bias_box(:,:,i), ...
        bias_gss(:,:,i), di(:,i), di_box(:,i), di_gss(:,i), ...
        di0(:,i), mag_bias_box(:,i), mag_bias_gss(:,i), ~, ~, ~] = ...
        impulse_resol(l, vr, u0, min_fres, max_fres, fres_inc, origin, props, ...
            err_level, num_ite, window_params, false);
end

% Average over data collected at the same feature resolution over different
% radii and show variation.
mean_dI0 = mean(dI0, 3);
sd_dI0 = std(dI0, 0, 3);
mean_bias_box = mean(bias_box, 3);
sd_bias_box = std(bias_box, 0, 3);
mean_bias_gss = mean(bias_gss, 3);
sd_bias_gss = std(bias_gss, 0, 3);

mean_mag_dI0 = mean(di0, 2);
mag_dI0_sd = std(di0, 0, 2);
mean_mag_bias_box = mean(mag_bias_box, 2);
mag_bias_box_sd = std(mag_bias_box, 0, 2);
mean_mag_bias_gss = mean(mag_bias_gss, 2);
mag_bias_gss_sd = std(mag_bias_gss, 0, 2);

mean_dI = mean(dI, 3);
dI_sd = std(dI, 0, 3);
mean_dI_box = mean(dI_box, 3);
dI_sd_box = std(dI_box, 0, 3);
mean_dI_gss = mean(dI_gss, 3);
dI_sd_gss = std(dI_gss, 0, 3);

mean_mag_dI = mean(di, 2);
mag_dI_sd = std(di, 0, 2);
mean_mag_dI_box = mean(di_box, 2);
mag_dI_sd_box = std(di_box, 0, 2);
mean_mag_dI_gss = mean(di_gss, 2);
mag_dI_sd_gss = std(di_gss, 0, 2);


% The plots mirror those in impulse_resol.m, only collecting from different
% radii.

%%%%%%%%%%%%%%%% Dimensional Plots %%%%%%%%%%%%%%%%%

% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [];
dim_str = {'x', 'y', 'z'};

% Font size for titles.
titleFsize = 11;

fres_ini = 1;
fresp = fres(fres_ini: end);

for dim = dims
    % Smoother bias plot.
    figure;
    errorbar(fresp, mean_dI0(dim, fres_ini:end), sd_dI0(dim, fres_ini:end), 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    errorbar(fresp, mean_bias_box(dim, fres_ini:end), sd_bias_box(dim, fres_ini:end), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    errorbar(fresp, mean_bias_gss(dim, fres_ini:end), sd_bias_gss(dim, fres_ini:end), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
    hold on

    legend({'imperfect resolution', ...
        'box filtered', ...
        'Gaussian-filtered'}, ...  
        'Interpreter', 'latex')
    xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
    ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
    title(strcat('$', dim_str{dim}, '$ Smoother Bias'), 'FontSize', titleFsize)

    % Mean error plot.
    figure;
    errorbar(fresp, mean_dI(dim, fres_ini:end), dI_sd(dim, fres_ini:end), 'ko', 'MarkerFaceColor','black', 'LineWidth', 1)
    hold on
    errorbar(fresp, mean_dI_box(dim, fres_ini:end), dI_sd_box(dim, fres_ini:end), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
    hold on
    errorbar(fresp, mean_dI_gss(dim, fres_ini:end), dI_sd_gss(dim, fres_ini:end), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)
    hold on

    legend({'unfiltered', ...
        'box-filtered', ...
        'Gaussian-filtered'}, ...  
        'Interpreter', 'latex')
    xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
    ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
    title(sprintf('$%s$ mean error at $\\delta u = %.0f\\%%$', string(dim_str{dim}), ...
        props(end)*100), 'FontSize', titleFsize)
end

%%%%%%%%%%%%%%%%%%% Magnitude Plots %%%%%%%%%%%%%%%%%%%%%

% Smoother bias plot.
figure;
errorbar(fresp, mean_mag_dI0(fres_ini:end), mag_dI0_sd(fres_ini:end), 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
hold on
errorbar(fresp, mean_mag_bias_box(fres_ini:end), mag_bias_box_sd(fres_ini:end), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
errorbar(fresp, mean_mag_bias_gss(fres_ini:end), mag_bias_gss_sd(fres_ini:end), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
hold on

legend({'imperfect resolution', ...
    'box filtered', ...
    'Gaussian-filtered'}, ...  
    'Interpreter', 'latex')
xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
ylabel('$\left|\frac{\delta I}{I}\right|$')
title('Resolution Errors', 'FontSize', titleFsize)

% Mean error plot.
figure;
errorbar(fresp, mean_mag_dI(fres_ini:end), mag_dI_sd(fres_ini:end), 'ko', 'MarkerFaceColor','black', 'LineWidth', 1)
hold on
errorbar(fresp, mean_mag_dI_box(fres_ini:end), mag_dI_sd_box(fres_ini:end), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
errorbar(fresp, mean_mag_dI_gss(fres_ini:end), mag_dI_sd_gss(fres_ini:end), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)
hold on

legend({'unfiltered', ...
    'box-filtered', ...
    'Gaussian-filtered'}, ...  
    'Interpreter', 'latex')
xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
ylabel('$\left|\frac{\delta I}{I}\right|$')
title(sprintf('Mean Error Magnitude at $\\delta u = %.0f$\\%%', ...
        props(end)*100), 'FontSize', titleFsize)
    
% Relative smoother bias plot to baseline resolution error.
figure;
scatter(fresp, abs(mean_mag_dI0(fres_ini:end) - mean_mag_bias_box(fres_ini:end)), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
scatter(fresp, abs(mean_mag_dI0(fres_ini:end) - mean_mag_bias_gss(fres_ini:end)), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
hold on

legend({'box filtered', ...
    'Gaussian-filtered'}, ...  
    'Interpreter', 'latex')
xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
ylabel('$\left|\frac{\delta I}{I} - \frac{\delta I_0}{I}\right|$')
title('Smoother bias relative to imperfection of resolution', 'FontSize', titleFsize)