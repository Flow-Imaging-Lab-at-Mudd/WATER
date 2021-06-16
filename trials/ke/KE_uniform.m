% Presume 'sp' is defined.
% sp = 0.1;
[x, y, z, u, v, w, Mag] = hill_vortex_3D(sp, 1, 1, 1);
vf = VelocityField.import_grid_separate(x,y,z,u,v,w);

% Minimal and maximal volume dimensions.
vol_range = [floor(1/4*vf.getDims())' floor(1/2*vf.getDims())'];
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
for i = 1: num_ite
    % Random range.
    range = randRange(vf.dims, vol_range);
    vf.setRange(range)
    
    % Run script for KE error samrpling.
    [dK(:,i), dK_box(:,i), dK_gss(:,i), bias_box(i), bias_gss(i)] = ...
        KE_uniform_err_run(vf, range);

%     pause
%     close all
end

% Global smoother error as a function of kinetic energy plot.
figure;
scatter(abs(dK(:)), abs(dK_box(:)), 'r', 'filled')
hold on
scatter(abs(dK(:)), abs(dK_gss(:)), 'b', 'filled')

legend('box-filtered', 'Gaussian-filtered', 'Interpreter', 'latex')
xlabel('Unfiltered $\left|\frac{\delta K}{K}\right|$')
ylabel('Filtered $\left|\frac{\delta K}{K}\right|$')
title(strcat('Normalized Spacing:', {' '}, string(sp)))
