% Uniform resolutions considered.
resol = flip(0.01: 0.05: 1);
resol_count = length(resol);

% Randomly sample effective regions.
num_ite = 20;

% Containers for data across all runs.
err_box = zeros(resol_count, num_ite);
err_gss = zeros(resol_count, num_ite);

bias_box = zeros(resol_count, num_ite);
bias_gss = zeros(resol_count, num_ite);

alpha_box = zeros(resol_count, num_ite);
alpha_gss = zeros(resol_count, num_ite);

% Introduce noise proportionally.
props = 0: 0.1: 1.5;

for k = 1: resol_count
    % Construct Hill vortex with specified resolution.
    [x, y, z, u, v, w, Mag] = hill_vortex_3D(resol(k), 1, 1, 1);
    vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
    vf.data.speed = sqrt(sum(vf.U.^2, 4));
    % Minimal and maximal volume dimensions.
    vol_range = [floor(1/4*vf.getDims())' 1/2*floor(vf.getDims())'];

    for i = 1: num_ite
        % This section mirrors KE_uniform.
        % Randomly subset range.
        range = randRange(vf.dims, vol_range);
        vf.setRange(range)
        % Run script for KE error samrpling.
        [err_box(k, i), err_gss(k, i), bias_box(k, i), ...
            bias_gss(k, i), alpha_box(k, i), alpha_gss(k, i)] = ...
            KE_uniform_resol_run(vf, range, props);
    end
end

mean_err_box = mean(err_box, 2);
mean_err_gss = mean(err_gss, 2);

smoother_bias_box = mean(bias_box, 2);
smoother_bias_gss = mean(bias_gss, 2);

amp_box = mean(alpha_box, 2);
amp_gss = mean(alpha_gss, 2);

figure;
errorbar(flip(resol), flip(mean_err_box), std(err_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth',1)
hold on
errorbar(flip(resol), flip(mean_err_gss), std(err_gss, 0, 2), 'ko', 'MarkerFaceColor','yellow', 'LineWidth', 1)

% Mean filter error plot.
legend('box-filtered', 'Gaussian-filtered')
xlabel('normalized spacing')
ylabel('$\bar{\left|\frac{\delta K}{K}\right|}$')
title('Mean Error Percentage of Smoother over Range of Noise')

% Smoother bias plot.
figure;
errorbar(flip(resol), flip(smoother_bias_box), std(bias_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
errorbar(flip(resol), flip(smoother_bias_gss), std(bias_gss, 0, 2), 'ko', 'MarkerFaceColor','yellow', 'LineWidth', 1)

legend('box-filtered', 'Gaussian-filtered')
xlabel('normalized spacing')
ylabel('$\left|\frac{\delta K}{K}\right|$')
title('Smoother Bias')

% Amplification coefficient plot.
figure;
errorbar(flip(resol), flip(amp_box), std(alpha_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
errorbar(flip(resol), flip(amp_gss), std(alpha_gss, 0, 2), 'ko', 'MarkerFaceColor','yellow', 'LineWidth', 1)

legend('box-filtered', 'Gaussian-filtered')
xlabel('normalized spacing')
ylabel('$\alpha$')
title('Quadratic Amplification Coefficient')
