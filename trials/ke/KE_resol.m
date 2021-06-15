% Uniform resolutions considered.
resol = flip(0.1: 0.1: 1);

% Randomly sample effective regions.
num_ite = 20;

% Containers for error measures with varying resolutions.
mean_err_box = zeros(1, length(resol));
mean_err_gss = zeros(1, length(resol));

smoother_bias_box = zeros(1, length(resol));
smoother_bias_gss = zeros(1, length(resol));

amp_box = zeros(1, length(resol));
amp_gss = zeros(1, length(resol));

% Containers for taking average across different volumes.
err_box = zeros(1, length(num_ite));
err_gss = zeros(1, length(num_ite));

bias_box = zeros(1, length(num_ite));
bias_gss = zeros(1, length(num_ite));

alpha_box = zeros(1, length(num_ite));
alpha_gss = zeros(1, length(num_ite));


for k = 1: length(resol)
    % Construct Hill vortex with specified resolution.
    [x, y, z, u, v, w, Mag] = hill_vortex_3D(resol(k), 1, 1, 1);
    vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
    vf.data.speed = sqrt(sum(vf.U.^2, 4));
    % Minimal and maximal volume dimensions.
    vol_range = [floor(1/4*vf.getDims())' floor(1/2*vf.getDims())'];

    for i = 1: num_ite
        % This section mirrors KE_uniform.
        % Randomly subset range.
        range = randRange(vf.dims, vol_range);
        vf.setRange(range)
        % Run script for KE error samrpling.
        [err_box(i), err_gss(i), bias_box(i), ...
            bias_gss(i), alpha_box(i), alpha_gss(i)] = KE_uniform_resol_run(vf, range);
    end
    
    mean_err_box(k) = mean(err_box);
    mean_err_gss(k) = mean(err_gss);
    
    smoother_bias_box(k) = mean(bias_box);
    smoother_bias_gss(k) = mean(bias_gss);

    amp_box(k) = mean(alpha_box);
    amp_gss(k) = mean(alpha_gss);
end

figure;
scatter(flip(resol), flip(mean_err_box), 'filled')
hold on
scatter(flip(resol), flip(mean_err_gss), 'filled')

legend('box-filtered', 'Gaussian-filtered')
xlabel('normalized spacing')
ylabel('$\bar{\left|\frac{\delta K}{K}\right|}$')
title('Mean Error Percentage of Smoother over Range of Noise')

figure;
scatter(flip(resol), flip(smoother_bias_box), 'filled')
hold on
scatter(flip(resol), flip(smoother_bias_gss), 'filled')

legend('box-filtered', 'Gaussian-filtered')
xlabel('normalized spacing')
ylabel('$\left|\frac{\delta K}{K}\right|$')
title('Smoother Bias')

figure;
scatter(flip(resol), flip(amp_box), 'filled')
hold on
scatter(flip(resol), flip(amp_gss), 'filled')

legend('box-filtered', 'Gaussian-filtered')
xlabel('normalized spacing')
ylabel('$\alpha$')
title('Quadratic Amplification Coefficient')