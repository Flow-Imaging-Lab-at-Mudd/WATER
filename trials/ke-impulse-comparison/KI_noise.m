% Comparison of noise profiles of impulse and KE computations.
%
% Derek Li, March 2022

clear
startup;

% Constant paremeters.
l = 1;
vr = 1;
r = l*vr;
spr = 0.05;
u0 = 1;
[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, u0, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

winsize = 8;
overlap = 0.5;
window_params = [];

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

[dI, dI_box, dI_gss, dI0, bias_box, bias_gss, dI_sd, dI_sd_box, dI_sd_gss, ...
    di, di_box, di_gss, di0, mag_bias_box, mag_bias_gss, di_sd, di_sd_box, di_sd_gss, vf] = ...
    impulse_err_run(vf, props, origin, I0, num_ite, window_params, false);

[dK, dK_box, dK_gss, dK0, Kbias_box, Kbias_gss, dK_sd, dK_sd_box, dK_sd_gss, vf] = ...
    KE_err_run(vf, props, K0, KEf, num_ite, window_params, false);

% Save plots.
savePlot = 0;
if savePlot
    img_fdr = sprintf('%s\\trials\\ke-impulse-comparison\\window_stats\\window-size\\o=%.2f\\dx=%.4f\\', ...
        rootFolder, overlap, spr);
    if ~isfolder(img_fdr)
        mkdir(img_fdr);
    end
end

% Computation without filter.
figure;
errorbar(props, dK, dK_sd, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
% Plot either the z-component error of impulse or the total error magnitude.
errorbar(props, dI(3,:), dI_sd(3,:), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
% errorbar(props, di_gss, di_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
legend({'kinetic energy', 'impulse'})
xlabel('Windowing size')
ylabel('Normalized error')
title(sprintf('Unfiltered error with noise with $o=%.2f$', overlap))

if savePlot
    try
        saveas(gcf, sprintf('%serr-%.2f.fig', img_fdr, props(2)))
        saveas(gcf, sprintf('%serr-%.2f.jpg', img_fdr, props(2)))
    catch
    end
end

% Computation with filter.
figure;
errorbar(props, dK_gss, dK_sd_gss, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
% Plot either the z-component error of impulse or the total error magnitude.
errorbar(props, dI_gss(3,:), dI_sd_gss(3,:), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
% errorbar(props, di_gss, di_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
legend({'kinetic energy', 'impulse'})
xlabel('Windowing size')
ylabel('Normalized error')
title(sprintf('Filtered error with noise with $o=%.2f$', overlap))

if savePlot
    try
        saveas(gcf, sprintf('%serr-%.2f.fig', img_fdr, props(2)))
        saveas(gcf, sprintf('%serr-%.2f.jpg', img_fdr, props(2)))
    catch
    end
end