% Vary the window size used in downsampling and present it effect on error.

% Original vortex used.
sp = 0.02;
fr = 1;
u0 = 1;
[x, y, z, u, v, w, Mag] = Hill_Vortex(sp, fr, 1, 1, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

% Origin chosen.
origin = [0 0 0]';

% Zoom in on vortical region.
vf.setRangePosition(fr * repmat([-1 1], 3, 1))

% Uniform windows in x, y, z.
windows = 2.^(1: 5);
windows_count = size(windows, 2);

% Constant overlap used.
overlap = 0.5;

% Proportional noise.
props = [0 1.5];
props_count = size(props, 2);

% Containers for error data at different window sizes.
dI = zeros(3, windows_count);
dI_sd = zeros(3, windows_count);
dId = zeros(3, windows_count);
dI_box = zeros(3, windows_count);
dI_sd_box = zeros(3, windows_count);
dI_gss = zeros(3, windows_count);
dI_sd_gss = zeros(3, windows_count);
bias_box = zeros(3, windows_count);
bias_gss = zeros(3, windows_count);

% Minimum feature resolution to reduce error beneath an error level.
err_level = 0.1;
min_bias_fres = zeros(3, windows_count);
min_err_fres = zeros(3, windows_count);
% Maximum freature size tried.
max_fres = 12;
min_fres = 1;
% Increment.
fres_inc = 1;

for i = 1: windows_count
    winsize = windows(i);
    [dId(:,i), dI(:,i), dI_sd(:,i), dI_box(:,i), dI_sd_box(:,i), dI_gss(:,i), ...
        dI_sd_gss(:,i), bias_box(:,i), bias_gss(:,i), ~] = ...
        impulse_err_window_stats(vf, props, origin, fr, u0, winsize, overlap);
    % Compute required resolution to compensate for 
    [min_bias_fres(:,i), min_err_fres(:,i)] = ...
        impulse_window_resol(fr, winsize, overlap, min_fres:fres_inc:max_fres, ...
            origin, err_level, props);
end

% 'dId' retains sign.
abs_dId = abs(dId);

% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [2 1 3];
dim_str = {'x', 'y', 'z'};

% Save plots.
img_fdr = strcat('C:\Users\derek\flow\trials\impulse\windowing\window-size\o=', ...
    string(overlap), '\');
mkdir(img_fdr);

for dim = dims
    % Plot of sampling bias.
    figure;
    scatter(windows, dId(dim, :), 'k', 'filled')
    xlabel('Window Size $w$')
    ylabel(strcat('$\frac{\delta I_', string(dim_str{dim}), '}{I}$'))
    title(strcat('$', dim_str{dim}, '$ Downsampling Bias'))
    saveas(gcf, strcat(img_fdr, 'db-', dim_str{dim}, '.jpg'));
    
    % Plot of bias.
    figure;
    scatter(windows, abs_dId(dim,:), 'k', 'filled')
    hold on
    scatter(windows, bias_box(dim,:), 'r', 'filled')
    hold on
    scatter(windows, bias_gss(dim,:), 'b', 'filled')

    bias_line_d = polyfit(windows, abs_dId(dim,:), 1);
    bias_line_box = polyfit(windows, bias_box(dim,:), 1);
    bias_line_gss = polyfit(windows, bias_gss(dim,:), 1);

    hold on
    bias_fit_d = polyplot(bias_line_d, windows, 'k');
    hold on
    bias_fit_box = polyplot(bias_line_box, windows, 'r');
    hold on
    bias_fit_gss = polyplot(bias_line_gss, windows, 'b');
    
    legend(strcat('downsampling bias $\delta$'), ...
        strcat('box bias'),  ...
        strcat('Gaussian bias'), ...
        strcat('downsampling fit $r^2 = $', {' '}, string(cor(bias_fit_d, abs_dId(dim,:)))), ...
        strcat('box fit $r^2 = $', {' '}, string(cor(bias_fit_box, bias_box(dim,:)))), ...
        strcat('Gaussian fit $r^2 = $', {' '}, string(cor(bias_fit_gss, bias_gss(dim,:)))), ...
        'Interpreter', 'latex')
    xlabel('Window Size $w$')
    ylabel(strcat('$\kappa_', dim_str{dim}, '$'))
    title(strcat('$', dim_str{dim}, '$ downsampled bias at overlap $o=$', ...
        {' '}, string(overlap)))
    saveas(gcf, strcat(img_fdr, 'bias-', dim_str{dim}, '.jpg'));
    
    % Plot of mean absolute errors.
    figure;
    errorbar(windows, dI(dim,:), dI_sd(dim,:), 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    errorbar(windows, dI_box(dim,:), dI_sd_box(dim,:), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    errorbar(windows, dI_gss(dim,:), dI_sd_gss(dim,:), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
    
    dId_line = polyfit(windows, abs_dId(dim,:), 1);
    dI_line_box = polyfit(windows, dI_box(dim,:), 1);
    dI_line_gss = polyfit(windows, dI_gss(dim,:), 1);

    hold on
    dId_fit = polyplot(dId_line, windows, 'k');
    hold on
    dI_fit_box = polyplot(dI_line_box, windows, 'r');
    hold on
    dI_fit_gss = polyplot(dI_line_gss, windows, 'b');

    legend(strcat('downsampling bias $\delta$'), ...
        strcat('box bias'),  ...
        strcat('Gaussian bias'), ...
        strcat('downsampling fit $r^2 = $', {' '}, string(cor(dId_fit, abs_dId(dim,:)))), ...
        strcat('box fit $r^2 = $', {' '}, string(cor(dI_fit_box, dI_box(dim,:)))), ...
        strcat('Gaussian fit $r^2 = $', {' '}, string(cor(dI_fit_gss, dI_gss(dim,:)))), ...
        'Interpreter', 'latex')
    
    xlabel('Window Size $w$')
    ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
    title(strcat('$', dim_str{dim}, '$ mean smoother error at overlap $o=$', {' '}, string(overlap)))
    saveas(gcf, strcat(img_fdr, 'err-', dim_str{dim}, '.jpg'));
end

%%%%%%%%%%% Magnitude plots %%%%%%%%%%%%%%

mag_dId = sqrt(sum(dId.^2, 1));
mag_bias_box = sqrt(sum(bias_box.^2, 1));
mag_bias_gss = sqrt(sum(bias_gss.^2, 1));

mag_dI = sqrt(sum(dI.^2, 1));
mag_dI_sd = sqrt(sum(dI_sd.^2, 1));
mag_dI_box = sqrt(sum(dI_box.^2, 1));
mag_dI_sd_box = sqrt(sum(dI_sd_box.^2, 1));
mag_dI_gss = sqrt(sum(dI_gss.^2, 1));
mag_dI_sd_gss = sqrt(sum(dI_sd_gss.^2, 1));

% Plot of sampling bias.
figure;
scatter(windows, mag_dId, 'k', 'filled')
xlabel('Window Size $w$')
ylabel('$|\frac{\delta I}{I}|$')
title('Downsampled Bias Magnitude')
saveas(gcf, strcat(img_fdr, 'db-mag.jpg'));

% Plot of bias.
figure;
scatter(windows, mag_dId, 'k', 'filled')
hold on
scatter(windows, mag_bias_box, 'r', 'filled')
hold on
scatter(windows, mag_bias_gss, 'b', 'filled')

bias_line_d = polyfit(windows, mag_dId, 1);
bias_line_box = polyfit(windows, mag_bias_box, 1);
bias_line_gss = polyfit(windows, mag_bias_gss, 1);

hold on
bias_fit_d = polyplot(bias_line_d, windows, 'k');
hold on
bias_fit_box = polyplot(bias_line_box, windows, 'r');
hold on
bias_fit_gss = polyplot(bias_line_gss, windows, 'b');

legend(strcat('downsampling bias $\delta$'), ...
    strcat('box bias'),  ...
    strcat('Gaussian bias'), ...
    strcat('downsampling fit $r^2 = $', {' '}, string(cor(bias_fit_d, mag_dId))), ...
    strcat('box fit $r^2 = $', {' '}, string(cor(bias_fit_box, mag_bias_box))), ...
    strcat('Gaussian fit $r^2 = $', {' '}, string(cor(bias_fit_gss, mag_bias_gss))), ...
    'Interpreter', 'latex')
xlabel('Window Size $w$')
ylabel('$\kappa$')
title(strcat('Downsampling bias magnitude at overlap $o=$', {' '}, string(overlap)))
saveas(gcf, strcat(img_fdr, 'bias-mag.jpg'));
    
% Plot of mean absolute errors.
figure;
errorbar(windows, mag_dI, mag_dI_sd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
hold on
errorbar(windows, mag_dI_box, mag_dI_sd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
errorbar(windows, mag_dI_gss, mag_dI_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)

dId_line = polyfit(windows, mag_dId, 1);
dI_line_box = polyfit(windows, mag_dI_box, 1);
dI_line_gss = polyfit(windows, mag_dI_gss, 1);

hold on
dId_fit = polyplot(dId_line, windows, 'k');
hold on
dI_fit_box = polyplot(dI_line_box, windows, 'r');
hold on
dI_fit_gss = polyplot(dI_line_gss, windows, 'b');

legend(strcat('downsampling bias $\delta$'), ...
    strcat('box bias'),  ...
    strcat('Gaussian bias'), ...
    strcat('downsampling fit $r^2 = $', {' '}, string(cor(dId_fit, mag_dId))), ...
    strcat('box fit $r^2 = $', {' '}, string(cor(dI_fit_box, mag_dI_box))), ...
    strcat('Gaussian fit $r^2 = $', {' '}, string(cor(dI_fit_gss, mag_dI_gss))), ...
    'Interpreter', 'latex')

xlabel('Window Size $w$')
ylabel('$|\frac{\delta I}{I}|$')
title(strcat('Mean smoother error magnitude at overlap $o=$', {' '}, string(overlap)))
saveas(gcf, strcat(img_fdr, 'err-mag.jpg'));

%%%%%% Make graphs for required minimum feature resolution.%%%%%%
% Bias graph.
figure;
scatter(windows, min_bias_fres(3, :), 'k', 'filled')
hold on
scatter(windows, min_bias_fres(1, :), 'r', 'filled')
hold on
scatter(windows, min_bias_fres(2, :), 'b', 'filled')

legend({'unfiltered', 'box-filtered', 'Gaussian-filtered'}, 'Interpreter', 'latex')
xlabel('Window Size $w$')
ylabel('$\frac{r}{s}$')
title(strcat('Minimum feature resolution for', ...
    {' '}, string(err_level*100), '\% smoother bias'))

% Mean error graph.
figure;
scatter(windows, min_err_fres(3, :), 'k', 'filled')
hold on
scatter(windows, min_err_fres(1, :), 'r', 'filled')
hold on
scatter(windows, min_err_fres(2, :), 'b', 'filled')

legend({'unfiltered', 'box-filtered', 'Gaussian-filtered'}, 'Interpreter', 'latex')
xlabel('Window Size $w$')
ylabel('$\frac{r}{s}$')
title(strcat('Minimum feature resolution for', ...
    {' '}, string(err_level*100), '\% mean error', ' given', {' '}, ...
    string(100*props(1)), '-', string(100*props(end)), '\% noise'))