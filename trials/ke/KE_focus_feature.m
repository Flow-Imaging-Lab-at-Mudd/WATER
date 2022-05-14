function KE_focus_feature(sp, fr, props)
% 'sp' is normalized spacing of Hill's Vortex.
% 'fr' is the feature radius of the vortex.

[x, y, z, u, v, w, ~] = hill_vortex_3D(sp, fr, 1, 1);
vf = VelocityField.import_grid_separate(x,y,z,u,v,w);

% Focus on the global region with the freestream velocity subtracted.
vf.addVelocity(-vf.U(1,1,1,:));

% Randomly sample effective regions.
num_ite = 10;

% Containers for data over iterations.
dK = zeros(length(props), num_ite);
dK_box = zeros(length(props), num_ite);
dK_gss = zeros(length(props), num_ite);
bias_box = zeros(1, num_ite);
bias_gss = zeros(1, num_ite);

for i = 1: num_ite
    % Run script for KE error sampling.
    [dK(:,i), dK_box(:,i), dK_gss(:,i), bias_box(i), bias_gss(i)] = ...
        KE_err_run(vf, props);

end

% Global smoother error as a function of kinetic energy plot.
K_noise = abs(dK(:));

figure;
scatter(K_noise, abs(dK_box(:)), 'r', 'filled')
hold on
scatter(K_noise, abs(dK_gss(:)), 'b', 'filled')
hold on

% 1-1 line.
plot(K_noise, K_noise, 'black')
hold on

% Smoother biases.
mean_bias_box = mean(abs(bias_box));
sd_bias_box = std(abs(bias_box));
yline(mean_bias_box, '-', 'Color', 'r')

hold on
mean_bias_gss = mean(abs(bias_gss));
sd_bias_gss = std(abs(bias_gss));
yline(mean_bias_gss, '-', 'Color', 'b')

legend('box-filtered', 'Gaussian-filtered', ...
    '$y=x$ identity line', ...
    strcat('box bias $\kappa = $', string(mean_bias_box), '$\pm$', string(sd_bias_box)), ...
    strcat('Gaussian bias $\kappa = $', string(mean_bias_gss), '$\pm$', string(sd_bias_gss)), ...
    'Interpreter', 'latex')
xlabel('Unfiltered $\left|\frac{\delta K}{K}\right|$')
ylabel('Filtered $\left|\frac{\delta K}{K}\right|$')
