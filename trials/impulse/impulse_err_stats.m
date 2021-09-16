function [dI_mean, dI_sd, dI_mean_box, dI_sd_box, dI_mean_gss, dI_sd_gss, ...
    bias_box, bias_gss, dI0] = impulse_err_stats(vf, props, origin, fr, u0)
% For synthetic data set, impulse_err_stats(vf, props, origin, fr, u0). The
% physical parameters 'fr', 'u0' are solely used to compute the theoretical
% vortical impulse of the Hill's vortex.
% For experimental data set,
% impulse_err_stats(vf, props, origin, I0), where 'origin_ref' is the
% origin from which the impulse taken as correct is computed.
% 
% Extract sample statistics so that the error is represented at a lower
% dimension.

% Complete error data.
% Synthetic.
[dI, dI_box, dI_gss, bias_box, bias_gss] = impulse_err_run(vf, props, origin, fr, u0);
% % Experimental.
% [dI, dI_box, dI_gss, bias_box, bias_gss] = impulse_err_run(vf, props, origin, origin_ref);

% Consider only absolute error.
abs_dI = abs(dI);
dI_mean = mean(abs_dI, 2);
dI_sd = std(abs_dI, 0, 2);

abs_dI_box = abs(dI_box);
dI_mean_box = mean(abs_dI_box, 2);
dI_sd_box = std(abs_dI_box, 0, 2);

abs_dI_gss = abs(dI_gss);
dI_mean_gss = mean(abs_dI_gss, 2);
dI_sd_gss = std(abs_dI_gss, 0, 2);

% Error due to imperfect resolution (perfect being infinite) and origin
% selection.
dI0 = abs_dI(:,1);

bias_box = abs(bias_box);
bias_gss = abs(bias_gss);

bias_mag_box = norm(bias_box);
bias_mag_gss = norm(bias_gss);