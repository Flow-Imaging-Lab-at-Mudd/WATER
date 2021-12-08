% Apply impulse_resol.m at different radii to generate graphs of error
% versus feature resolution, varied both by global resolution and feature
% size.

% Constant freestream velocity.
u0 = 1;

% Define global spacing with feature resolution, the minimum of which is
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
mag_dI0 = zeros(sps_count, radii_count);

bias_box = zeros(3, sps_count, radii_count);
bias_gss = zeros(3, sps_count, radii_count);
mag_bias_box = zeros(sps_count, radii_count);
mag_bias_gss = zeros(sps_count, radii_count);


mag_dI = zeros(sps_count, radii_count);
mag_dI_box = zeros(sps_count, radii_count);
mag_dI_gss = zeros(sps_count, radii_count);

% Proportions of velocity noise introduced.
props = [0 2];

% Collect data for different radii.
for i = 1: radii_count
    fr = radii(i);
    % Vary global resolution for the given radius.
    [fres, dI(:,:,i), dI_box(:,:,i), dI_gss(:,:,i), bias_box(:,:,i), ...
        bias_gss(:,:,i), mag_dI(:,i), mag_dI_box(:,i), mag_dI_gss(:,i), ...
        mag_bias_box(:,i), mag_bias_gss(:,i), dI0(:,:,i), mag_dI0(:,i)] = ...
        impulse_resol(fr, u0, min_fres, max_fres, fres_inc, origin, props);
end

% Average over data collected at the same feature resolution over different
% radii and show variation.
mean_dI0 = mean(dI0, 3);
sd_dI0 = std(dI0, 0, 3);
mean_bias_box = mean(bias_box, 3);
sd_bias_box = std(bias_box, 0, 3);
mean_bias_gss = mean(bias_gss, 3);
sd_bias_gss = std(bias_gss, 0, 3);

mean_mag_dI0 = mean(mag_dI0, 2);
mag_dI0_sd = std(mag_dI0, 0, 2);
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

mean_mag_dI = mean(mag_dI, 2);
mag_dI_sd = std(mag_dI, 0, 2);
mean_mag_dI_box = mean(mag_dI_box, 2);
mag_dI_sd_box = std(mag_dI_box, 0, 2);
mean_mag_dI_gss = mean(mag_dI_gss, 2);
mag_dI_sd_gss = std(mag_dI_gss, 0, 2);


% The plots mirror those in impulse_resol.m, only collecting from different
% radii.

%%%%%%%%%%%%%%%% Dimensional Plots %%%%%%%%%%%%%%%%%

% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [2 1 3];
dim_str = {'x', 'y', 'z'};

% Font size for titles.
titleFsize = 11;

fres_ini = 2;
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
    title(strcat('$', string(dim_str{dim}), '$ Mean Error at $\delta u = \,$', ...
        string(props(end)*100), '\%'), 'FontSize', titleFsize)
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
title(strcat('Mean Error Magnitude at $\delta u = \,$', ...
        string(props(end)*100), '\%'), 'FontSize', titleFsize)
    
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