function [dK, dK_box, dK_gss, dK0, bias_box, bias_gss, dK_sd, dK_sd_box, dK_sd_gss, vf] = ...
    KE_err_run_constN(vf, props, K0, KEf, num_ite, window_params, display_plots)

if ~exist('display_plots', 'var')
    display_plots = 0;
end

[dK, dK_box, dK_gss, dK0, bias_box, bias_gss, dK_sd, dK_sd_box, dK_sd_gss, vf] = ...
    KE_err_run(vf, props, K0, KEf, num_ite, window_params, display_plots);
% Select first noise level as the only one considered.
dK = dK(2);
dK_box = dK_box(2);
dK_gss = dK_gss(2);
dK_sd = dK_sd(2);
dK_sd_box = dK_sd_box(2);
dK_sd_gss = dK_sd_gss(2);