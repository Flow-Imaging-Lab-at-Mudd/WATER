% Presume 'sp' (spacing), 'fr' (feature radius) are defined.
% sp = 0.05;
% fr = 1;
[x, y, z, u, v, w, Mag] = hill_vortex_3D(sp, fr, 1, 1);
vf = VelocityField.import_grid_separate(x,y,z,u,v,w);

% Span of feature in indices.
fs = 2 * floor(fr ./ vf.resol) + 1;
fs = min([fs; vf.getDims()], [], 1);
% Feature is at the center of the grid. Random range intersects the feature.
% Proportions of random region to central vortex feature.
region_props = [2/3 3/2];
vol_range = [floor(region_props(1)*fs)' ...
    floor(min([region_props(2)*fs; vf.getDims()-1], [], 1))'];
% Randomly sample effective regions.
num_ite = 10;
% Introduce noise proportionally.
props = 0: 0.1: 3;

% Containers for data over iterations.
dK = zeros(length(props), num_ite);
dK_box = zeros(length(props), num_ite);
dK_gss = zeros(length(props), num_ite);
bias_box = zeros(1, num_ite);
bias_gss = zeros(1, num_ite);

% vf.setRange([floor(vf.getDims()/2 - fs/2); floor(vf.getDims()/2 + fs/2)]')
% vf.plotVector(vf.U_e, 0, '')

for i = 1: num_ite
    % Random range.
    range = randRange(vf.dims, vol_range);
    span = range(:,2) - range(:,1);
    % Anchor with randomness. No longer vectorized here.
    shift_prop = 1/3;
    left_idx = floor(vf.getDims()/2 - fs/2) + randi(int32(1/2*shift_prop*[-fs(1) fs(1)]));
    % Ensure not exceeding index range.
    left_idx = max([ones(1, 3); left_idx])';
    right_idx = min([vf.getDims()' left_idx+span], [], 2);
    range = [left_idx right_idx];
    
%     % Ensure intersection.
%     vf.plotVector(vf.U_e, 0, '$\vec{u}$');
    
    % Run script for KE error sampling.
    [dK(:,i), dK_box(:,i), dK_gss(:,i), bias_box(i), bias_gss(i)] = ...
        KE_err_run(vf, range, props);

%     pause
%     close all
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
