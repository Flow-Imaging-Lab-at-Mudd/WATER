% Original vortex used.
sp = 0.02;
fr = 1;
[x, y, z, u, v, w, Mag] = hill_vortex_3D(sp, fr, 1, 1);
vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
% Subtract freestream velocity.
vf.addVelocity(-vf.U(1,1,1,:))

% Zoom in on vortical region.
vf.setRangePosition(fr * repmat([-1 1], 3, 1))

% Uniform windows in x, y, z.
windows = 2.^(1: 5);
windows_count = size(windows, 2);

% Constant overlap used.
overlap = 0.75;

% Proportional noise.
props = 0: 0.1: 3;
props_count = size(props, 2);

% Containers for error data at different window sizes.
dK = zeros(props_count, windows_count);
dKd = zeros(1, windows_count);
dK_box = zeros(props_count, windows_count);
dK_gss = zeros(props_count, windows_count);
bias_box = zeros(1, windows_count);
bias_gss = zeros(1, windows_count);
mean_abs_err_box = zeros(1, windows_count);
mean_abs_err_gss = zeros(1, windows_count);

% Minimum feature resolution to reduce error beneath an error level.
err_level = 0.1;
min_bias_fres = zeros(2, windows_count);
min_err_fres = zeros(2, windows_count);
% Feature resolutions tried.
max_fres = 12;
min_fres = 1;
% Increment.
fres_inc = 1;

for i = 1: windows_count
    winsize = windows(i);
    [dK(:,i), dKd(i), dK_box(:,i), dK_gss(:,i), bias_box(i), bias_gss(i), ...
        mean_abs_err_box(i), mean_abs_err_gss(i), ~] = ...
        KE_err_window_run(vf, winsize, overlap, props);
    % Compute the minimum feature resolution required.
    [min_bias_fres(:,i), min_err_fres(:,i)] = ...
        KE_window_resol(winsize, overlap, min_fres:fres_inc:max_fres, ...
        err_level, props);
end

% Plot of sampling bias, expected monotonic increase of underestimation.
figure;
scatter(windows, dKd, 'k', 'filled')
xlabel('Window Size $w$')
ylabel('$\kappa$')
title('Downsampling Bias')

% Plot of bias.
figure;
scatter(windows, abs(dKd), 'k', 'filled')
hold on
scatter(windows, abs(bias_box), 'r', 'filled')
hold on
scatter(windows, abs(bias_gss), 'b', 'filled')

bias_line_d = polyfit(windows, abs(dKd), 1);
bias_line_box = polyfit(windows, abs(bias_box), 1);
bias_line_gss = polyfit(windows, abs(bias_gss), 1);

hold on
bias_fit_d = polyplot(bias_line_d, windows, 'k');
hold on
bias_fit_box = polyplot(bias_line_box, windows, 'r');
hold on
bias_fit_gss = polyplot(bias_line_gss, windows, 'b');


legend(strcat('downsampling bias $\delta$'), ...
    strcat('box bias'),  ...
    strcat('Gaussian bias'), ...
    strcat('downsampling fit $r^2 = $', {' '}, string(cor(bias_fit_d, abs(dKd)))), ...
    strcat('box fit $r^2 = $', {' '}, string(cor(bias_fit_box, abs(bias_box)))), ...
    strcat('Gaussian fit $r^2 = $', {' '}, string(cor(bias_fit_gss, abs(bias_gss)))), ...
    'Interpreter', 'latex')
xlabel('Window Size $w$')
ylabel('$\kappa$')
title(strcat('Downsampling biases at overlap $o=$', {' '}, string(overlap)))

% Save plots.
img_fdr = strcat('C:\Users\derek\flow\trials\ke\windowing\window-size\o=', ...
    string(overlap), '\');
mkdir(img_fdr);
% saveas(gcf, strcat(img_fdr, 'bias-', string(overlap), 'o.jpg'));

% Plot of mean absolute errors.
figure;
scatter(windows, mean_abs_err_box, 'r', 'filled')
hold on
scatter(windows, mean_abs_err_gss, 'b', 'filled')

legend('box-filtered', 'Gaussian-filtered')
xlabel('Window Size $w$')
ylabel('$\left|\frac{\delta K}{K}\right|$')
title(strcat('Mean Smoother Errors at Overlap $o=$', {' '}, string(overlap)))
% saveas(gcf, strcat(img_fdr, 'err-', string(overlap), 'o.jpg'));

%%%%%% Make graphs for required minimum feature resolution.
% Bias graph.
figure;
scatter(windows, min_bias_fres(1, :), 'r', 'filled')
hold on
scatter(windows, min_bias_fres(2, :), 'b', 'filled')

legend({'box-filtered', 'Gaussian-filtered'}, 'Interpreter', 'latex')
xlabel('Window Size $w$')
ylabel('$\frac{r}{s}$')
title(strcat('Minimum feature resolution for', ...
    {' '}, string(err_level*100), '\% smoother bias'))

% Mean error graph.
figure;
scatter(windows, min_err_fres(1, :), 'r', 'filled')
hold on
scatter(windows, min_err_fres(2, :), 'b', 'filled')

legend({'box-filtered', 'Gaussian-filtered'}, 'Interpreter', 'latex')
xlabel('Window Size $w$')
ylabel('$\frac{r}{s}$')
title(strcat('Minimum feature resolution for', ...
    {' '}, string(err_level*100), '\% mean error', ' given', {' '}, ...
    string(100*props(1)), '-', string(100*props(end)), '\% noise'))