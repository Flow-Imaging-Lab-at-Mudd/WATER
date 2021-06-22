% Original vortex used.
sp = 0.02;
fr = 1;
[x, y, z, u, v, w, Mag] = hill_vortex_3D(sp, fr, 1, 1);
vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
% Subtract freestream velocity.
vf.addVelocity(-vf.U(1,1,1,:))

% Uniform windows in x, y, z.
windows = 2.^(1: 6);
windows_count = size(windows, 2);

% Constant overlap used.
overlap = 0.75;

% Proportional noise.
props = 0.1: 0.1: 3;
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

for i = 1: windows_count
    winsize = windows(i);
    [dK(:,i), dKd(i), dK_box(:,i), dK_gss(:,i), bias_box(i), bias_gss(i), ...
        mean_abs_err_box(i), mean_abs_err_gss(i)] = ...
        KE_err_window_run(vf, winsize, overlap, props, fr);
end

% Plot of sampling bias, expected monotonic increase of underestimation.
figure;
scatter(windows, dKd, 'k', 'filled')

% Plot of bias.
figure;
scatter(windows, abs(dKd), 'k', 'filled')
hold on
scatter(windows, abs(bias_box), 'r', 'filled')
hold on
scatter(windows, abs(bias_gss), 'b', 'filled')

legend(strcat('downsampling bias $\delta$'), ...
    strcat('box bias'),  ...
    strcat('Gaussian bias'), 'Interpreter', 'latex')
xlabel('Window Size $w$')
ylabel('$\kappa$')
title(strcat('Downsampling biases at overlap $o=$', {' '}, string(overlap)))

% Save plots.
img_fdr = strcat('C:\Users\derek\flow\trials\ke\windowing\window-size\o=', ...
    string(overlap), '\');
mkdir(img_fdr);
saveas(gcf, strcat(img_fdr, 'bias-', string(overlap), 'o.jpg'));

% Plot of mean absolute errors.
figure;
scatter(windows, mean_abs_err_box, 'r', 'filled')
hold on
scatter(windows, mean_abs_err_gss, 'b', 'filled')

legend('box-filtered', 'Gaussian-filtered')
xlabel('Window Size $w$')
ylabel('$\left|\frac{\delta K}{K}\right|$')
title(strcat('Mean Smoother Errors at Overlap $o=$', {' '}, string(overlap)))
saveas(gcf, strcat(img_fdr, 'err-', string(overlap), 'o.jpg'));