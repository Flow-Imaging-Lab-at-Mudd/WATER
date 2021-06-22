% Original vortex used.
sp = 0.05;
fr = 1;
[x, y, z, u, v, w, Mag] = hill_vortex_3D(sp, fr, 1, 1);
vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
% Subtract freestream velocity.
vf.addVelocity(-vf.U(1,1,1,:))

% Uniform windows.
windows = 2: 11;
windows_count = size(windows, 2);

% Constant overlap used.
overlap = 0.75;

% Proportional noise.
props = 0.1: 0.1: 3;
props_count = size(props, 2);

% Containers for error data at different window sizes.
dK = zeros(props_count, windows_count);
dKd = zeros(props_count, windows_count);
bias = zeros(1, windows_count);
mean_err = zeros(1, windows_count);

for i = 1: windows_count
    winsize = windows(i);
    [dK(:,i), dKd(:,i), bias(i), mean_err(i)] = ...
        KE_err_window_run(vf, winsize, overlap, props, fr);
end

% Plot of bias.
figure;
scatter(windows, abs(bias), 'r', 'filled')
xlabel('Window Size $w$')
ylabel('$\kappa$')
title(strcat('Downsampling biases at overlap $o=$', {' '}, string(overlap)))

% Save plots.
img_fdr = strcat('C:\Users\derek\flow\trials\ke\windowing\window-size\o=', ...
    string(overlap), '\');
mkdir(img_fdr);
saveas(gcf, strcat(img_fdr, 'bias-', string(overlap), 'o.jpg'));

% Plot of mean errors.
figure;
scatter(windows, mean_err, 'r', 'filled')
xlabel('Window Size $w$')
ylabel('$\left|\frac{\delta K}{K}\right|$')
title(strcat('Mean smoother errors at overlap $o=$', {' '}, string(overlap)))
saveas(gcf, strcat(img_fdr, 'err-', string(overlap), 'o.jpg'));