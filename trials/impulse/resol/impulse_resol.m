function [fres, dI, dI_box, dI_gss, bias_box, bias_gss, mag_dI, mag_dI_box, ...
    mag_dI_gss, mag_bias_box, mag_bias_gss, dI0, mag_dI0] = ...
    impulse_resol(fr, u0, min_fres, max_fres, fres_inc, origin, props, show_plot)
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

% Generate range of spacings for evenly spaced feature
% resolutions.
% Global spacing of downsampled data.
sps = zeros(1, floor((max_fres-1)/fres_inc) + 1);
% Minimal feature resolution.
sps(1) = fr / min_fres;

for i = 2: size(sps, 2)
    sps(i) = fr / (fres_inc + fr/sps(i-1));
end

sps_count = size(sps, 2);

% Feature resolution.
fres = fr ./ sps;

% Consider stochastic effect.
num_ite = 10;

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
    sp = sps(k);
    [x, y, z, u, v, w, ~] = Hill_Vortex(sp, fr, u0, 1, 1);
    vf = VelocityField.importCmps(x, y, z, u, v, w);
    
    % Focus on vortical region.
    vf.setRangePosition(fr*repmat([-1 1], 3, 1))
    
    for i = 1: num_ite
        % Run script for impulse error sampling.
        [err(:,k,i), err_sd(:,k,i), err_box(:,k,i), err_sd_box(:,k,i), ...
            err_gss(:,k,i), err_sd_gss(:,k,i), bias_box(:,k), ...
            bias_gss(:,k), dI0(:,k)] = ...
            impulse_err_stats(vf, props, origin, ...
                HillImpulse(vf.fluid.density, vf.scale.len, fr, u0));
    end
end

% Average over trials.
dI = mean(err, 3);
dI_sd = mean(err_sd, 3);
dI_box = mean(err_box, 3);
dI_sd_box = mean(err_sd_box, 3);
dI_gss = mean(err_gss, 3);
dI_sd_gss = mean(err_sd_gss, 3);

mag_dI0 =  sqrt(sum(dI0.^2, 1));
mag_bias_box = sqrt(sum(bias_box.^2, 1));
mag_bias_gss = sqrt(sum(bias_gss.^2, 1));

% Magnitudes.
mag_dI = sqrt(sum(dI.^2, 1));
mag_dI_sd = sqrt(sum(dI_sd.^2, 1));
mag_dI_box = sqrt(sum(dI_box.^2, 1));
mag_dI_sd_box = sqrt(sum(dI_sd_box.^2, 1));
mag_dI_gss = sqrt(sum(dI_gss.^2, 1));
mag_dI_sd_gss = sqrt(sum(dI_sd_gss.^2, 1));


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
        'box filtered', ...
        'Gaussian-filtered'}, ...
        'Interpreter', 'latex')
    xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
    ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
    title(strcat('$', dim_str{dim}, '$ Smoother Bias at $r = $', {' '}, string(fr)))
    
    % Mean error plot.
    figure;
    errorbar(fres, dI(dim,:), dI_sd(dim,:), 'ko', 'MarkerFaceColor','black', 'LineWidth', 1)
    hold on
    errorbar(fres, dI_box(dim,:), dI_sd_box(dim,:), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
    hold on
    errorbar(fres, dI_gss(dim,:), dI_sd_gss(dim,:), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)
    
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
scatter(fres, mag_dI0, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
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
errorbar(fres, mag_dI, mag_dI_sd, 'ko', 'MarkerFaceColor','black', 'LineWidth', 1)
hold on
errorbar(fres, mag_dI_box, mag_dI_sd_box, 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
errorbar(fres, mag_dI_gss, mag_dI_sd_gss, 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)
hold on

legend({'unfiltered', ...
    'box-filtered', ...
    'Gaussian-filtered'}, ...  
    'Interpreter', 'latex')
xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
ylabel('$\left|\frac{\delta I}{I}\right|$')
title(strcat('Mean Error Magnitude at $\delta u = $', ...
        string(props(2)*100), '\% at $r = $', {' '}, string(fr)))