function [dI_mean, dI_sd, dI_mean_box, dI_sd_box, dI_mean_gss, dI_sd_gss, ...
    bias_box, bias_gss, dI0] = impulse_err_stats(vf, props, origin, I0)
% See impulse_err_run.m for a description of the computation performed. This
% function is a wrapper which extracts the average of absolute errors over
% the different levels of noises, if applicable.
%
% Derek Li, March 2022

% Complete error data.
[dI, dI_box, dI_gss, bias_box, bias_gss] = impulse_err_run(vf, props, origin, I0, 1);

% Consider only absolute error.
abs_dI = abs(dI);
% Average only over noisy results.
dI_mean = mean(abs_dI(:, 2:end), 2);
dI_sd = std(abs_dI(:, 2:end), 0, 2);

abs_dI_box = abs(dI_box);
dI_mean_box = mean(abs_dI_box(:, 2:end), 2);
dI_sd_box = std(abs_dI_box(:, 2:end), 0, 2);

abs_dI_gss = abs(dI_gss);
dI_mean_gss = mean(abs_dI_gss(:, 2:end), 2);
dI_sd_gss = std(abs_dI_gss(:, 2:end), 0, 2);

% Error due to imperfect resolution (perfect being infinite) and origin
% selection.
dI0 = abs_dI(:,1);

bias_box = abs(bias_box);
bias_gss = abs(bias_gss);