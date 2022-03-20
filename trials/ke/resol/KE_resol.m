function [fres, dK, dK_box, dK_gss, bias_box, bias_gss, dK0, dK_sd, dK_box_sd, dK_gss_sd] = ...
    KE_resol(l, vr, u0, min_fres, max_fres, fres_inc, props, err_level, display_plots)

% Radius of vortex.
fr = l*vr;
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
% Normalize spacing.
spr = sps / fr;

num_ite = 10;
% Containers for data across all runs.
% Errors here are mean absolute errors.
% Errors here are mean absolute errors.
err = zeros(sps_count, num_ite);
err_box = zeros(sps_count, num_ite);
err_gss = zeros(sps_count, num_ite);
% Standard deviations are also the average over trials.
err_sd = zeros(sps_count, num_ite);
err_sd_box = zeros(sps_count, num_ite);
err_sd_gss = zeros(sps_count, num_ite);

% Resolution errors.
bias_box = zeros(sps_count, 1);
bias_gss = zeros(sps_count, 1);
dK0 = zeros(sps_count, 1);

% Compute vortical KE of Hill's vortex.
KEf = @(vf, n) Hill_VortKE(vf, fr, n);
% Vortex parameters.
density = 1000;
len_unit = 1e-3;
% Theoretical KE.
K0 = Hill_KE(density, len_unit, fr, u0);

for k = 1: sps_count
    % Construct Hill vortex with specified resolution.
    [x, y, z, u, v, w] = Hill_Vortex(spr(k), l, vr, u0, 1);
    vf = VelocityField.importCmps(x, y, z, u, v, w);
    % Focus on vortical region.
    vf.setRangePosition(fr*repmat([-1 1], 3, 1))
    
    % Compute error per resolution.
    for i = 1: num_ite
        % Run script for impulse error sampling.
        [err(k,i), ~, err_box(k,i), ~, err_gss(k,i), ~, ...
            bias_box(k), bias_gss(k), dK0(k)] = KE_err_stats(vf, props, K0, KEf);
    end
end

% Average over trials.
dK = mean(err, 2);
dK_sd = std(err, 0, 2);
dK_box = mean(err_box, 2);
dK_box_sd = std(err_box, 0, 2);
dK_gss = mean(err_gss, 2);
dK_gss_sd = std(err_gss, 0, 2);

%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%
if ~exist('display_plots', 'var') || ~display_plots
    return
end

% Smoother bias plot.
figure;
scatter(fres, dK0, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
hold on
scatter(fres, bias_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
scatter(fres, bias_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)

legend({'imperfect resolution', 'box-filtered', 'Gaussian-filtered'})
xlabel(strcat('Feature Resolution'))
ylabel(strcat('$\left|\frac{\delta K}{K}\right|$'))
title(strcat('Resolution errors at $r = $', {' '}, string(fr)))

% Mean error plot.
figure;
errorbar(fres, dK, dK_sd, 'ko', 'MarkerFaceColor','black', 'LineWidth', 1)
hold on
errorbar(fres, dK_box, dK_box_sd, 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
errorbar(fres, dK_gss, dK_gss_sd, 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)

legend({'unfiltered', ...
    'box-filtered', ...
    'Gaussian-filtered'})
xlabel(strcat('Feature Resolution'))
ylabel(strcat('$\left|\frac{\delta K}{K}\right|$'))
title(sprintf('Mean Error at $\\delta u = %.0f$\\%% at $r = %.0f$', string(props(2)*100), string(fr)))

% Compute minimum resolution needed when smoothers are applied. First row
% corresponds to biases; second row, noisy error. Box, first column;
% Gaussian, second.
min_res = -1*ones(2, 2);

% Lowest feature resolutions to achieve the desired error level.
try
    min_res(1,1) = fres(end - find(flip(bias_box) < err_level, 1) + 1);
catch
end
try
    min_res(1,2) = fres(end - find(flip(bias_gss) < err_level, 1) + 1);
catch
end

try
    min_res(2,1) = fres(end - find(flip(dK_box) < err_level, 1) + 1);
catch
end
try
    min_res(2,2) = fres(end - find(flip(dK_gss) < err_level, 1) + 1);
catch
end

% Print resolution requirements.
fprintf('Minimum resolution for %.0f %% bias: \n', err_level*100)
fprintf('Box: %d \n', min_res(1,1))
fprintf('Gaussian: %d \n', min_res(1,2))

fprintf('Minimum resolution for %.0f %% noise-propagated error: \n', err_level*100)
fprintf('Box: %d \n', min_res(2,1))
fprintf('Gaussian: %d \n', min_res(2,2))

