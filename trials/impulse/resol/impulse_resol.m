function [fres, dI, dI_box, dI_gss, bias_box, bias_gss, di, di_box, ...
    di_gss, mag_bias_box, mag_bias_gss, dI0, di0, di_sd, di_box_sd, di_gss_sd] = ...
    impulse_resol(l, vr, u0, min_fres, max_fres, fres_inc, origin, props, show_plot)
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

% Global spacing of downsampled data.
sps = zeros(1, floor((max_fres-1)/fres_inc) + 1);
% Minimal feature resolution.
sps(1) = fr / min_fres;
% Generate range of spacings for evenly spaced feature
% resolutions.
for i = 2: size(sps, 2)
    sps(i) = fr / (fres_inc + fr/sps(i-1));
end

sps_count = size(sps, 2);
% Feature resolution.
fres = fr ./ sps;
% Normalize spacing.
spr = sps / fr;

% Consider stochastic effect.
num_ite = 10;

% Vortex parameters.
density = 1000;
len_unit = 1e-3;
% Theoretical impulse values.
I0 = Hill_Impulse(density, len_unit, fr, u0);

% Containers for data across all runs.
% Errors here are mean absolute errors.
err = zeros(3, sps_count, num_ite);
err_box = zeros(3, sps_count, num_ite);
err_gss = zeros(3, sps_count, num_ite);
% Standard deviations are also the average over trials.
err_sd = zeros(3, sps_count, num_ite);
err_sd_box = zeros(3, sps_count, num_ite);
err_sd_gss = zeros(3, sps_count, num_ite);

bias_box = zeros(3, sps_count, 1);
bias_gss = zeros(3, sps_count, 1);

% Error due purely to the imperfection of resolution and origin selection.
dI0 = zeros(3, sps_count);

for k = 1: sps_count
    % Construct Hill vortex with specified resolution.
    [x, y, z, u, v, w] = Hill_Vortex(spr(k), l, vr, u0, 1);
    vf = VelocityField.importCmps(x, y, z, u, v, w);
    
    % Focus on vortical region.
    vf.setRangePosition(fr*repmat([-1 1], 3, 1))
    
    for i = 1: num_ite
        % Run script for impulse error sampling.
        [err(:,k,i), ~, err_box(:,k,i), ~, err_gss(:,k,i), ~, bias_box(:,k), ...
            bias_gss(:,k), dI0(:,k)] = impulse_err_stats(vf, props, origin, I0);
    end
end

% Magnitudes.
di = squeeze(sqrt(sum(err.^2, 1)));
di_box = squeeze(sqrt(sum(err_box.^2, 1)));
di_gss = squeeze(sqrt(sum(err_gss.^2, 1)));

% Average over trials.
di_sd = std(di, 0, 2);
di_box_sd = std(di_box, 0, 2);
di_gss_sd = std(di_gss, 0, 2);
dI_sd = std(err, 0, 3);
dI_box_sd = std(err_box, 0, 3);
dI_gss_sd = std(err_gss, 0, 3);

dI = mean(err, 3);
dI_box = mean(err_box, 3);
dI_gss = mean(err_gss, 3);
di = mean(di, 2);
di_box = mean(di_box, 2);
di_gss = mean(di_gss, 2);

di0 = sqrt(sum(dI0.^2, 1));
mag_bias_box = sqrt(sum(bias_box.^2, 1));
mag_bias_gss = sqrt(sum(bias_gss.^2, 1));

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