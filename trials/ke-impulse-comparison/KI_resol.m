% Comparison of resolution profiles of impulse and KE computations.

clear
startup;

% Constant parameters.
l = 1;
vr = 1;
r = l*vr;
u0 = 1;

% Windowing parameters.
winsize = 16;
overlap = 0.5;
window_params = [winsize overlap];

% Noise levels.
props = [0 1.5];
% Desired level of error.
err_level = 0.1;
% Iterations.
num_ite = 5;

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

% Store the VFs with proper resolutions generated in the KE trials.
[~, dK, dK_box, dK_gss, Kbias_box, Kbias_gss, dK0, dK_sd, dK_box_sd, dK_gss_sd, vfds] = ...
    KE_resol(l, vr, u0, min_fres, max_fres, fres_inc, props, err_level, num_ite, window_params, false);

[fres, dI, dI_box, dI_gss, dI0, Ibias_box, Ibias_gss, di, di_box, ...
    di_gss, di0, mag_bias_box, mag_bias_gss, dI_sd, dI_box_sd, dI_gss_sd,...
    di_sd, di_box_sd, di_gss_sd] = ...
    impulse_resol(l, vr, u0, min_fres, max_fres, fres_inc, origin, props, err_level, num_ite, window_params, false, vfds);

% Resolutions plotted.
fres_ini = 2;
fresp = fres(fres_ini: end);

% Save plots.
savePlot = 0;
if savePlot
    img_fdr = sprintf('%s\\trials\\ke-impulse-comparison\\window_stats\\resolutions\\w=%d\\o=%.2f\\', ...
        rootFolder, winsize, overlap);
    if ~isfolder(img_fdr)
        mkdir(img_fdr);
    end
end

% Baseline resolution error.
figure;
scatter(fresp, dK0(fres_ini:end), 'filled', 'r')
hold on
% Plot either the z-component error of impulse or the total error magnitude.
% scatter(fresp, di0(fres_ini:end), 'filled', 'b')
scatter(fresp, dI0(3,fres_ini:end), 'filled', 'b')
legend({'kinetic energy', 'impulse'})
xlabel('Feature resolution')
ylabel('Normalized error')
title('Baseline resolution error')

if savePlot
    saveas(gcf, sprintf('%sres_err.fig', img_fdr))
    saveas(gcf, sprintf('%sres_err.jpg', img_fdr))
end

% Resolution error by filter.
figure;
scatter(fresp, Kbias_gss(fres_ini:end), 'filled', 'r')
hold on
% Plot either the z-component error of impulse or the total error magnitude.
% scatter(fresp, mag_bias_gss(fres_ini:end), 'filled', 'b')
scatter(fresp, Ibias_gss(3,fres_ini:end), 'filled', 'b')
legend({'kinetic energy', 'impulse'})
xlabel('Feature resolution')
ylabel('Normalized error')
title('Resolution error of filter')

if savePlot
    saveas(gcf, sprintf('%sfilter_err.fig', img_fdr))
    saveas(gcf, sprintf('%sfilter_err.jpg', img_fdr))
end


% Error with noise.
figure;
errorbar(fresp, dK_gss(fres_ini:end), dK_gss_sd(fres_ini:end), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
% Plot either the z-component error of impulse or the total error magnitude.
% errorbar(fresp, di_gss(fres_ini:end), di_gss_sd(fres_ini:end), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
errorbar(fresp, dI_gss(3,fres_ini:end), dI_gss_sd(3,fres_ini:end), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
legend({'kinetic energy', 'impulse'})
xlabel('Feature resolution')
ylabel('Normalized error')
title('Resolution error with noise')

if savePlot
    try
        saveas(gcf, sprintf('%serr-%.2f.fig', img_fdr, props(2)))
        saveas(gcf, sprintf('%serr-%.2f.jpg', img_fdr, props(2)))
    catch
    end
end
