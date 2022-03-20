% Comparison of noise profiles of impulse and KE computations.
%
% Derek Li, March 2022

% Constant paremeters.
l = 1;
vr = 1;
r = l*vr;
spr = 0.05;
u0 = 1;
[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, u0, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

% Noise levels.
props = 0:0.2:3;
% Iterations.
num_ite = 10;

% Theoretical values.
origin = [0 0 0]';
I0 = Hill_Impulse(vf.fluid.density, vf.scale.len, r, u0, r);
% Compute vortical KE of Hill's vortex.
KEf = @(vf, n) Hill_VortKE(vf, r, n);
% Theoretical KE.
K0 = Hill_KE(vf.fluid.density, vf.scale.len, r, u0);

[dI, dI_box, dI_gss, Ib_box, Ib_gss, did, did_box, did_gss] = ...
    impulse_err_run(vf, props, origin, I0, num_ite, false);
% Consider error magnitude for comparison.
di = squeeze(sqrt(sum(dI.^2, 1)));
di_box = squeeze(sqrt(sum(dI_box.^2, 1)));
di_gss = squeeze(sqrt(sum(dI_gss.^2, 1)));

[dK, dK_box, dK_gss, Kb_box, Kb_gss, dKd, dKd_box, dKd_gss] = ...
    KE_err_run(vf, props, K0, KEf, num_ite, false);
dK = abs(dK);
dK_box = abs(dK_box);
dK_gss = abs(dK_gss);
Kb_box = abs(Kb_box);
Kb_gss = abs(Kb_gss);


% Unfiltered error.
figure;
errorbar(props, dK, dKd, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
errorbar(props, di, did, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
legend({'kinetic energy', 'impulse'})
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('Normalized error')
title('Unfiltered errors of KE vs. impulse computations')

% Box-filtered error.
figure;
errorbar(props, dK_box, dKd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
errorbar(props, di_box, did_box, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
legend({'kinetic energy', 'impulse'})
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('Normalized error')
title('Box-filtered errors of KE vs. impulse computations')

% Unfiltered error.
figure;
errorbar(props, dK_gss, dKd_gss, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
errorbar(props, di_gss, did_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
legend({'kinetic energy', 'impulse'})
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('Normalized error')
title('Gaussian-filtered errors of KE vs. impulse computations')
