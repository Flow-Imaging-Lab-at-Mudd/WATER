% Original vortex used.
sp = 0.02;
fr = 1;
[x, y, z, u, v, w, Mag] = hill_vortex_3D(sp, fr, 1, 1);
vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
% Subtract freestream velocity.
vf.addVelocity(-vf.U(1,1,1,:))

% Constant window size in x, y, z.
winsize = 16;

% Constant overlap used.
overlap = 0.25: 0.25: 0.75;
overlap_count = size(overlap, 2);

% Proportional noise.
props = 0.1: 0.1: 3;
props_count = size(props, 2);

% Containers for error data at different window sizes.
dK = zeros(props_count, overlap_count);
dKd = zeros(1, overlap_count);
dK_box = zeros(props_count, overlap_count);
dK_gss = zeros(props_count, overlap_count);
bias_box = zeros(1, overlap_count);
bias_gss = zeros(1, overlap_count);
mean_abs_err_box = zeros(1, overlap_count);
mean_abs_err_gss = zeros(1, overlap_count);

for i = 1: overlap_count
    op = overlap(i);
    [dK(:,i), dKd(i), dK_box(:,i), dK_gss(:,i), bias_box(i), bias_gss(i), ...
        mean_abs_err_box(i), mean_abs_err_gss(i), ~] = ...
        KE_err_window_run(vf, winsize, op, props, fr);
end

% Plot of bias.
figure;
scatter(overlap, abs(dKd), 'k', 'filled')
hold on
scatter(overlap, abs(bias_box), 'r', 'filled')
hold on
scatter(overlap, abs(bias_gss), 'b', 'filled')

bias_line_d = polyfit(overlap, abs(dKd), 1);
bias_line_box = polyfit(overlap, abs(bias_box), 1);
bias_line_gss = polyfit(overlap, abs(bias_gss), 1);

hold on
bias_fit_d = polyplot(bias_line_d, overlap, 'k');
hold on
bias_fit_box = polyplot(bias_line_box, overlap, 'r');
hold on
bias_fit_gss = polyplot(bias_line_gss, overlap, 'b');


legend(strcat('downsampling bias $\delta$'), ...
    strcat('box bias'),  ...
    strcat('Gaussian bias'), ...
    strcat('downsampling fit $r^2 = $', {' '}, string(cor(bias_fit_d, abs(dKd)))), ...
    strcat('box fit $r^2 = $', {' '}, string(cor(bias_fit_box, abs(bias_box)))), ...
    strcat('Gaussian fit $r^2 = $', {' '}, string(cor(bias_fit_gss, abs(bias_gss)))), ...
    'Interpreter', 'latex')
xlabel('Proportion of Window Overlap $o$')
ylabel('$\kappa$')
title(strcat('Downsampling Biases at $w=$', {' '}, string(winsize)))

% Save plots.
img_fdr = strcat('C:\Users\derek\flow\trials\ke\windowing\overlap\w=', ...
    string(winsize), '\');
mkdir(img_fdr);
saveas(gcf, strcat(img_fdr, 'bias-', string(winsize), 'w.jpg'));

% Plot of mean absolute errors.
figure;
scatter(overlap, mean_abs_err_box, 'r', 'filled')
hold on
scatter(overlap, mean_abs_err_gss, 'b', 'filled')

legend('box-filtered', 'Gaussian-filtered')
xlabel('Proportion of Window Overlap $o$')
ylabel('$\left|\frac{\delta K}{K}\right|$')
title(strcat('Mean Smoother Errors at $w=$', {' '}, string(winsize)))
saveas(gcf, strcat(img_fdr, 'err-', string(winsize), 'w.jpg'));