% Comparison of resolution profiles of impulse and KE computations.
%
% Derek Li, March 2022

% Constant parameters.
l = 1;
vr = 1;
r = l*vr;
u0 = 1;

% Noise levels.
props = [0 1.5];
% Desired level of error.
err_level = 0.1;
% Iterations.
num_ite = 10;

% Resolutions.
min_fres = 1;
fres_inc = 1;
max_fres = 10;

% Vortex parameters.
density = 1000;
len_unit = 1e-3;

% Theoretical values.
origin = [0 0 0]';
I0 = Hill_Impulse(density, len_unit, r, u0, r);
% Compute vortical KE of Hill's vortex.
KEf = @(vf, n) Hill_VortKE(vf, r, n);
% Theoretical KE.
K0 = Hill_KE(density, len_unit, r, u0);

[~, dK, dK_box, dK_gss, bias_box, bias_gss, dK0, dK_sd, dK_box_sd, dK_gss_sd] = ...
    KE_resol(l, vr, u0, min_fres, max_fres, fres_inc, props, err_level, false);

[fres, ~, ~, ~, ~, ~, di, di_box, ...
    di_gss, mag_bias_box, mag_bias_gss, dI0, di0, di_sd, di_box_sd, di_gss_sd] = ...
    impulse_resol(l, vr, u0, min_fres, max_fres, fres_inc, origin, props, false);

% Resolutions plotted.
fres_ini = 2;
fresp = fres(fres_ini: end);

% Resolution errors.
figure;
scatter(fresp, bias_gss(fres_ini:end), 'filled', 'r')
hold on
scatter(fresp, mag_bias_gss(fres_ini:end), 'filled', 'b')
legend({'kinetic energy', 'impulse'})
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('Normalized error')
title('Baseline resolution errors')

% Errors with noise.
figure;
errorbar(fresp, dK_gss(fres_ini:end), dK_gss_sd(fres_ini:end), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
errorbar(fresp, di_gss(fres_ini:end), di_gss_sd(fres_ini:end), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
legend({'kinetic energy', 'impulse'})
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('Normalized error')
title('Resolution errors with noise')
