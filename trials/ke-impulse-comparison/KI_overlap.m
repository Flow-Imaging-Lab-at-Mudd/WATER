% Comparison of error profiles of impulse and KE computations with respect
% to windowing overlap.


clear
startup;

% Constant parameters.
l = 1;
vr = 1;
r = l*vr;
u0 = 1;
spr = 0.025;

[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, u0, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

% Constant parameters.
fr = l*vr;
% Theoretical values of impulse & KE.
I0 = Hill_Impulse(vf.fluid.density, vf.scale.len, fr, u0);
K0 = Hill_KE(vf.fluid.density, vf.scale.len, fr, u0);
% Origin used for impulse calculation.
origin = [0 0 0]';
% Compute vortival KE of Hill's vortex.
KEf = @(vf, n) Hill_VortKE(vf, 1, n);

% Windowing parameters.
winsize = 16;
overlaps = 0.25*(1:3);
% Show the dimensions of downsampled field, assuming the field has uniform
% dimensions.
downsampled_dim(vf.dims(1), winsize, overlaps);

% Noise level.
props = [0 1.5];

% Compute KE profile.
[dK, dK_box, dK_gss, dK0, Kbias_box, Kbias_gss, dK_sd, dK_sd_box, dK_sd_gss] = ...
    KE_overlap(vf, K0, KEf, props, winsize, overlaps);

% Compute impulse profile.
[dI, dI_box, dI_gss, dI0, Ibias_box, Ibias_gss, ...
    di, di_box, di_gss, di0, ibias_box, ibias_gss, ...
    dI_sd, dI_sd_box, dI_sd_gss, di_sd, di_sd_box, di_sd_gss] = ...
    impulse_overlap(vf, I0, origin, props, winsize, overlaps);

% Save plots.
savePlot = 1;
if savePlot
    img_fdr = sprintf('%s\\trials\\ke-impulse-comparison\\window_stats\\overlap\\w=%d\\dx=%.4f\\', ...
        rootFolder, winsize, spr);
    if ~isfolder(img_fdr)
        mkdir(img_fdr);
    end
end

% Baseline resolution error.
figure;
scatter(overlaps, dK0, 'filled', 'r')
hold on
% Plot either the z-component error of impulse or the total error magnitude.
scatter(overlaps, dI0(3,:), 'filled', 'b')
% scatter(overlaps, di0, 'filled', 'b')
xticks(overlaps)
legend({'kinetic energy', 'impulse'})
xlabel('Overlap ratio of window')
ylabel('Normalized error')
title(sprintf('Baseline windowing error with $w=%d$', winsize))

if savePlot
    saveas(gcf, sprintf('%sres_err.fig', img_fdr))
    saveas(gcf, sprintf('%sres_err.jpg', img_fdr))
end

% Resolution error by filter.
figure;
scatter(overlaps, Kbias_gss, 'filled', 'r')
hold on
% Plot either the z-component error of impulse or the total error magnitude.
scatter(overlaps, Ibias_gss(3,:), 'filled', 'b')
% scatter(overlaps, ibias_gss, 'filled', 'b')
xticks(overlaps)
legend({'kinetic energy', 'impulse'})
xlabel('Overlap ratio of window')
ylabel('Normalized error')
title(sprintf('Baseline error of filter with $w=%d$', winsize))

if savePlot
    saveas(gcf, sprintf('%sfilter_err.fig', img_fdr))
    saveas(gcf, sprintf('%sfilter_err.jpg', img_fdr))
end

% Error with noise.
figure;
errorbar(overlaps, dK_gss, dK_sd_gss, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
% Plot either the z-component error of impulse or the total error magnitude.
errorbar(overlaps, dI_gss(3,:), dI_sd_gss(3,:), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
% errorbar(overlaps, di_gss, di_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
xticks(overlaps)
legend({'kinetic energy', 'impulse'})
xlabel('Overlap ratio of window')
ylabel('Normalized error')
title(sprintf('Error with noise with $w=%d$', winsize))

if savePlot
    try
        saveas(gcf, sprintf('%serr-%.2f.fig', img_fdr, props(2)))
        saveas(gcf, sprintf('%serr-%.2f.jpg', img_fdr, props(2)))
    catch
    end
end
