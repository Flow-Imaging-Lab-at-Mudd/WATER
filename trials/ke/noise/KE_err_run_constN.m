function [dK, dK_box, dK_gss, dK0, bias_box, bias_gss, dK_sd, dK_sd_box, dK_sd_gss, vf, ...
    dKm, dKm_box, dKm_gss, dKm_sd, dKm_sd_box, dKm_sd_gss] = ...
    KE_err_run_constN(vf, props, K0, KEf, num_ite, window_params, display_plots)
% A wrapper for 'KE_err_run.m' for the special case that only one noise
% level is specified, i.e. props = [0 n], where n is the level of noise in
% proportion to the mean speed. The values returned by this function are
% scalars corresponding to this noise level.

if ~exist('display_plots', 'var')
    display_plots = {};
end

[dK, dK_box, dK_gss, dK0, bias_box, bias_gss, dK_sd, dK_sd_box, dK_sd_gss, vf, ...
    dKm, dKm_box, dKm_gss, dKm_sd, dKm_sd_box, dKm_sd_gss] = ...
    KE_err_run(vf, props, K0, KEf, num_ite, window_params, display_plots);
% Select first noise level (second entry, the first being for no noise) as the only one considered.
dK = dK(2);
dK_box = dK_box(2);
dK_gss = dK_gss(2);
dK_sd = dK_sd(2);
dK_sd_box = dK_sd_box(2);
dK_sd_gss = dK_sd_gss(2);

dKm = dKm(2);
dKm_box = dKm_box(2);
dKm_gss = dKm_gss(2);
dKm_sd = dKm_sd(2);
dKm_sd_box = dKm_sd_box(2);
dKm_sd_gss = dKm_sd_gss(2);