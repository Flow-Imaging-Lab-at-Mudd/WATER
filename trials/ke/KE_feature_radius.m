% Uniform resolutions considered.
sp = 0.02;
radii = flip(0.05: 0.05: 1);
radii_count = length(radii);

% Randomly sample effective regions.
num_ite = 1;

% Containers for data across all runs.
err_box = zeros(radii_count, num_ite);
err_gss = zeros(radii_count, num_ite);

bias_box = zeros(radii_count, num_ite);
bias_gss = zeros(radii_count, num_ite);

alpha_box = zeros(radii_count, num_ite);
alpha_gss = zeros(radii_count, num_ite);

% Introduce noise proportionally.
props = 0: 0.1: 1.5;

for k = 1: radii_count
    % Construct Hill vortex with specified resolution.
    fr = radii(k);
    [x, y, z, u, v, w, Mag] = hill_vortex_3D(sp, fr, 1, 1);
    vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
    vf.data.speed = sqrt(sum(vf.U.^2, 4));
    % Span of feature in indices.
    fs = 2 * floor(fr ./ vf.resol) + 1;
    fs = min([fs; vf.getDims()], [], 1);
%     Feature is at the center of the grid. Random range intersects the feature.
%     Proportions of random region to central vortex feature.
    region_props = [2/3 3/2];
    vol_range = [floor(region_props(1)*fs)' ...
        floor(min([region_props(2)*fs; vf.getDims()-1], [], 1))'];

    for i = 1: num_ite
        % This section mirrors KE_focus_feature.
        % Random range.
        range = randRange(vf.dims, vol_range);
        span = range(:,2) - range(:,1);
        % Anchor with randomness. No longer vectorized here.
        shift_prop = 1/3;
        left_idx = floor(vf.getDims()/2 - fs/2) + randi(int32(1/2*shift_prop*[-fs(1) fs(1)]));
        left_idx = floor(vf.getDims()/2 - fs/2);
        % Ensure not exceeding index range.
        left_idx = max([ones(1, 3); left_idx])';
        right_idx = min([vf.getDims()' left_idx+floor(fs/2)], [], 2);
%         disp('--------')
%         disp(['radius'])
%         radii(k)
        range = [floor(vf.getDims()/2 - floor(fs/2)); floor(vf.getDims()/2 + floor(fs/2))]';
        vf.setRange([left_idx right_idx])
        
        % Run script for KE error samrpling.
        [err_box(k, i), err_gss(k, i), bias_box(k, i), ...
            bias_gss(k, i), alpha_box(k, i), alpha_gss(k, i)] = ...
            KE_mean_err_run(vf, props);
    end
end

mean_err_box = mean(err_box, 2);
mean_err_gss = mean(err_gss, 2);

smoother_bias_box = mean(bias_box, 2);
smoother_bias_gss = mean(bias_gss, 2);

amp_box = mean(alpha_box, 2);
amp_gss = mean(alpha_gss, 2);

figure;
errorbar(flip(radii), flip(mean_err_box), std(err_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth',1)
hold on
errorbar(flip(radii), flip(mean_err_gss), std(err_gss, 0, 2), 'ko', 'MarkerFaceColor','yellow', 'LineWidth', 1)

% Mean filter error plot.
legend('box-filtered', 'Gaussian-filtered')
xlabel('nomalized feature radius')
ylabel('$\bar{\left|\frac{\delta K}{K}\right|}$')
title('Mean Error Percentage of Smoother over Range of Noise')

% Smoother bias plot.
figure;
errorbar(flip(radii), flip(smoother_bias_box), std(bias_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
errorbar(flip(radii), flip(smoother_bias_gss), std(bias_gss, 0, 2), 'ko', 'MarkerFaceColor','yellow', 'LineWidth', 1)

legend('box-filtered', 'Gaussian-filtered')
xlabel('nomalized feature radius')
ylabel('$\left|\frac{\delta K}{K}\right|$')
title('Smoother Bias')

% Amplification coefficient plot.
figure;
errorbar(flip(radii), flip(amp_box), std(alpha_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
errorbar(flip(radii), flip(amp_gss), std(alpha_gss, 0, 2), 'ko', 'MarkerFaceColor','yellow', 'LineWidth', 1)

legend('box-filtered', 'Gaussian-filtered')
xlabel('nomalized feature radius')
ylabel('$\alpha$')
title('Quadratic Amplification Coefficient')
