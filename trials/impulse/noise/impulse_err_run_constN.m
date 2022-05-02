function [dI, dI_box, dI_gss, dI0, bias_box, bias_gss, dI_sd, dI_sd_box, dI_sd_gss, ...
    di, di_box, di_gss, di0, mag_bias_box, mag_bias_gss, di_sd, di_sd_box, di_sd_gss, vf] = ...
    impulse_err_run_constN(vf, props, origin, I0, num_ite, window_params, display_plots)
% A wrapper for 'impulse_err_run.m' for the special case that only one noise
% level is specified, i.e. props = [0 n], where n is the level of noise in
% proportion to the mean speed.

if ~exist('display_plots', 'var')
    display_plots = 0;
end

[dI, dI_box, dI_gss, dI0, bias_box, bias_gss, dI_sd, dI_sd_box, dI_sd_gss, ...
    di, di_box, di_gss, di0, mag_bias_box, mag_bias_gss, di_sd, di_sd_box, di_sd_gss, vf] = ...
    impulse_err_run(vf, props, origin, I0, num_ite, window_params, display_plots);

% Select first noise level as the only one considered.
dI = dI(:,2);
dI_box = dI_box(:,2);
dI_gss = dI_gss(:,2);
dI_sd = dI_sd(:,2);
dI_sd_box = dI_sd_box(:,2);
dI_sd_gss = dI_sd_gss(:,2);

di = di(2);
di_box = di_box(2);
di_gss = di_gss(2);
di_sd = di_sd(2);
di_sd_box = di_sd_box(2);
di_sd_gss = di_sd_gss(2);