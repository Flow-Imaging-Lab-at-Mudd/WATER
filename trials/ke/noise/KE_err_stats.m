function [dK_mean, dK_sd, dK_mean_box, dK_sd_box, dK_mean_gss, dK_sd_gss, ...
    bias_box, bias_gss, dK0] = KE_err_stats(vf, props, K0, KEf)
% See KE_err_run.m for a description of function. This function is a
% wrapper which extracts the average of absolute errors over the different
% levels of noises, if applicable.
%
% Derek Li, March 2022

% All errors are given in absolute normalized value.
[dK, dK_box, dK_gss, bias_box, bias_gss] = KE_err_run(vf, props, K0, KEf, 1);

abs_dK = abs(dK);
% Average only over noisy results.
dK_mean = mean(abs_dK(:, 2:end), 2);
dK_sd = std(abs_dK(:, 2:end), 0, 2);

abs_dK_box = abs(dK_box);
dK_mean_box = mean(abs_dK_box(:, 2:end), 2);
dK_sd_box = std(abs_dK_box(:, 2:end), 0, 2);

abs_dK_gss = abs(dK_gss);
dK_mean_gss = mean(abs_dK_gss(:, 2:end), 2);
dK_sd_gss = std(abs_dK_gss(:, 2:end), 0, 2);

% Error due to imperfect resolution.
dK0 = abs_dK(:,1);

bias_box = abs(bias_box);
bias_gss = abs(bias_gss);