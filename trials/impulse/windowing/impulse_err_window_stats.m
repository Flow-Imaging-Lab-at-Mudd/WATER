function [dId, dI_mean, dI_sd, dI_mean_box, dI_sd_box, dI_mean_gss, ...
    dI_sd_gss, bias_box, bias_gss, vfd] = ...
    impulse_err_window_stats(vf0, props, origin, fr, u0, winsize, overlap)
% Identical to impulse_err_window_run.m except that the means and standard
% deviations of absolute erros are extracted.

[dI, dId, dI_box, dI_gss, bias_box, bias_gss, vfd] = ...
    impulse_err_window_run(vf0, props, origin, fr, u0, winsize, overlap);

% Consider only absolute error.
abs_dI = abs(dI);
dI_mean = mean(abs_dI(2:end), 2);
dI_sd = std(abs_dI(2:end), 0, 2);

abs_dI_box = abs(dI_box);
dI_mean_box = mean(abs_dI_box(2:end), 2);
dI_sd_box = std(abs_dI_box(2:end), 0, 2);

abs_dI_gss = abs(dI_gss);
dI_mean_gss = mean(abs_dI_gss(2:end), 2);
dI_sd_gss = std(abs_dI_gss(2:end), 0, 2);

bias_box = abs(bias_box);
bias_gss = abs(bias_gss);

bias_mag_box = norm(bias_box);
bias_mag_gss = norm(bias_gss);