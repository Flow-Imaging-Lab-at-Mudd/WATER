% Uniform resolutions considered.
resol = flip(0.01: 0.05: 1);
resol_count = length(resol);

% Constant feature radius used while varying global resolution.
fr = 0.7;

% Randomly sample effective regions.
num_ite = 10;

% Containers for data across all runs.
err_box = zeros(resol_count, num_ite);
err_gss = zeros(resol_count, num_ite);

bias_box = zeros(resol_count, 1);
bias_gss = zeros(resol_count, 1);

alpha_box = zeros(resol_count, num_ite);
alpha_gss = zeros(resol_count, num_ite);

% Introduce noise proportionally.
props = 0: 0.1: 1.5;

for k = 1: resol_count
    % Construct Hill vortex with specified resolution.
    [x, y, z, u, v, w, Mag] = hill_vortex_3D(resol(k), fr, 1, 1);
    vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
    
    % Subtract freestream velocity to focus on central feature region.
    vf.addVelocity(-vf.U(1,1,1,:))
    
    for i = 1: num_ite
        % This section mirrors KE_uniform.
        % Run script for KE error samrpling.
        [err_box(k, i), err_gss(k, i), bias_box(k), ...
            bias_gss(k), alpha_box(k, i), alpha_gss(k, i)] = ...
            KE_mean_err_run(vf, props, fr);
    end
end

mean_err_box = mean(err_box, 2);
mean_err_gss = mean(err_gss, 2);

smoother_bias_box = bias_box;
smoother_bias_gss = bias_gss;

amp_box = mean(alpha_box, 2);
amp_gss = mean(alpha_gss, 2);

figure;
errorbar(flip(resol), flip(mean_err_box), std(err_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth',1)
hold on
errorbar(flip(resol), flip(mean_err_gss), std(err_gss, 0, 2), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)

% Mean filter error plot.
legend('box-filtered', 'Gaussian-filtered')
xlabel('Normalized Spacing $s$')
ylabel('$\bar{\left|\frac{\delta K}{K}\right|}$')
title(strcat('Mean Error of Smoother over $\delta u = $', {' '}, ...
    string(props(1)), '-', string(props(end)*100), '\% $\bar{u}$', ...
    ' at $r = $', {' '}, string(fr)))

% Folder for plots.
img_fdr = strcat('C:\Users\derek\flow\trials\ke\resol\');
mkdir(strcat(img_fdr, 'mean-error\'));
mkdir(strcat(img_fdr, 'smoother-bias'));
mkdir(strcat(img_fdr, 'amp-coeff'))
saveas(gcf, strcat(img_fdr, '\mean-error\err-', string(fr), 'r', '.jpg'));

% Smoother bias plot.
figure;
scatter(flip(resol), flip(smoother_bias_box), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
scatter(flip(resol), flip(smoother_bias_gss), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)

legend('box-filtered', 'Gaussian-filtered')
xlabel('Normalized Spacing $s$')
ylabel('$\left|\frac{\delta K}{K}\right|$')
title(strcat('Smoother Bias at $r = $', {' '}, string(fr)))

saveas(gcf, strcat(img_fdr, 'smoother-bias\bias-', string(fr), 'r', '.jpg'));

% Amplification coefficient plot.
figure;
errorbar(flip(resol), flip(amp_box), std(alpha_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
errorbar(flip(resol), flip(amp_gss), std(alpha_gss, 0, 2), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)

legend('box-filtered', 'Gaussian-filtered')
xlabel('Normalized Spacing $s$')
ylabel('$\alpha$')
title('Quadratic Amplification Coefficient')

saveas(gcf, strcat(img_fdr, 'amp-coeff\amp-', string(fr), 'r', '.jpg'));
