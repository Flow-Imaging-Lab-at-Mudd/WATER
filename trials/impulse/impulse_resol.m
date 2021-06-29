function [fres, dI, dI_box, dI_gss, mag_dI, mag_dI_box, mag_dI_gss] = ...
    impulse_resol(fr, max_fres, fres_inc, origin, props)

% % Constant feature radius used while varying global resolution.
% fr = 0.4;
% % Reference point for impulse calculation.
% origin = [0 0 0]';
% 
% props = 0: 0.1: 3;
% 
% max_fres = 20;
% fres_inc = 1;

% Custom freestream velocity.
u0 = 1;

% Generate range of spacings for evenly spaced feature
% resolutions.
% Global spacing of downsampled data.
sps = zeros(1, floor((max_fres-1)/fres_inc) + 1);
% Minimal feature resolution of 1.
sps(1) = fr;

for i = 2: size(sps, 2)
    sps(i) = fr / (fres_inc + fr/sps(i-1));
end

% Reverse ordering of spacing so that high resolutions precedes low.
sps = flip(sps);
sps_count = size(sps, 2);

% Feature resolution.
fres = fr ./ sps;

% Consider stochastic effect.
num_ite = 5;

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

% Flag for under-resolution.
under_resolved = false;

for k = 1: sps_count
    % Construct Hill vortex with specified resolution.
    sp = sps(k);
    [x, y, z, u, v, w, ~] = hill_vortex_3D(sp, fr, u0, 1);
    vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
    
    % Subtract freestream velocity to focus on central feature region.
    vf.addVelocity(-vf.U(1,1,1,:))
    % Focus on vortical region.
    vf.setRangePosition(fr*repmat([-1 1], 3, 1))
    
    
    for i = 1: num_ite
        % Run script for KE error samrpling.
        [err(:,k,i), err_sd(:,k,i), err_box(:,k,i), err_sd_box(:,k,i), ...
            err_gss(:,k,i), err_sd_gss(:,k,i), bias_box(:,k), bias_gss(:,k)] = ...
            impulse_err_stats(vf, props, origin, fr, u0);
    end
end

% Average over trials.
dI = mean(err, 3);
dI_sd = mean(err_sd, 3);
dI_box = mean(err_box, 3);
dI_sd_box = mean(err_sd_box, 3);
dI_gss = mean(err_gss, 3);
dI_sd_gss = mean(err_sd_gss, 3);


%%%%%%%%%%%%%%%% Dimensional Plots %%%%%%%%%%%%%%%%%
% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [2 1 3];
dim_str = {'x', 'y', 'z'};

% for dim = dims
%     % Smoother bias plot.
%     figure;
%     scatter(fres, bias_box(dim,:), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
%     hold on
%     scatter(fres, bias_gss(dim,:), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)
%     hold on
% 
%     legend({'box filtered', ...
%         'Gaussian-filtered'}, ...  
%         'Interpreter', 'latex')
%     xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
%     ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
%     title(strcat('$', dim_str{dim}, '$ Smoother Bias at $r = $', {' '}, string(fr)))
% 
%     % Mean error plot.
%     figure;
%     errorbar(fres, dI(dim,:), dI_sd(dim,:), 'ko', 'MarkerFaceColor','black', 'LineWidth', 1)
%     hold on
%     errorbar(fres, dI_box(dim,:), dI_sd_box(dim,:), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
%     hold on
%     errorbar(fres, dI_gss(dim,:), dI_sd_gss(dim,:), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)
%     hold on
% 
%     legend({'unfiltered', ...
%         'box-filtered', ...
%         'Gaussian-filtered'}, ...  
%         'Interpreter', 'latex')
%     xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
%     ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
%     title(strcat('$', string(dim_str{dim}), '$ Mean Error over $\delta u = $', ...
%         string(props(1)*100), '-', string(props(end)*100), '\% at $r = $', {' '}, string(fr)))
% end

%%%%%%%%%%%%%%%%%%% Magnitude Plots %%%%%%%%%%%%%%%%%%%%%
mag_bias_box = sqrt(sum(bias_box.^2, 1));
mag_bias_gss = sqrt(sum(bias_gss.^2, 1));

mag_dI = sqrt(sum(dI.^2, 1));
mag_dI_sd = sqrt(sum(dI_sd.^2, 1));
mag_dI_box = sqrt(sum(dI_box.^2, 1));
mag_dI_sd_box = sqrt(sum(dI_sd_box.^2, 1));
mag_dI_gss = sqrt(sum(dI_gss.^2, 1));
mag_dI_sd_gss = sqrt(sum(dI_sd_gss.^2, 1));

% Smoother bias plot.
figure;
scatter(fres, mag_bias_box, 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
scatter(fres, mag_bias_gss, 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)
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
title(strcat('Mean Error Magnitude over $\delta u = $', ...
        string(props(1)*100), '-', string(props(end)*100), '\% at $r = $', {' '}, string(fr)))