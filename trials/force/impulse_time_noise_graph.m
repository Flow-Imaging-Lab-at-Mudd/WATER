% Load data from folder.
fdr = 'C:\Users\derek\flow\data\L18_Run1\';
data = dir(fdr);

% Subfolder names corresponding to frames.
frame_fdrs = strings(1, size(data, 1));
% Import subfolders in sorted order.
for i = 1: size(data, 1)
    frame_fdrs(i) = strcat(data(i).name, '\');
end
frame_fdrs = sort(frame_fdrs);
% First two entries not valid subfolders. Next two with null velocity
% fields, so with last two.
frame_fdrs = frame_fdrs(5:end-2);

% Number of time frames to consider. Maximal 36.
num_frame = 36;
% Array of velocity fields, corresonding to different time frames.
vfs = cell(1, num_frame);

% Time progression.
t = 1/2000*(1: num_frame);

% Fixed origin, center of vortex in the first frame.
center = [-8 10 -21]';
% Fixed position range where impulse is calculated.
vort = [-30 10; center(2)-15 center(2)+15; -35 -5];

% Impulse computed in said region.
I = zeros(3, num_frame);

% Load data as velocity fields.
for i = 1: num_frame
    % Post-processed data.
    load(strcat(fdr, frame_fdrs(i), '\3DPIV_postprocessed_results_calibrated.mat'), ...
        'x', 'y', 'z', 'u', 'v', 'w')
    vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
    % Focus on region.
    vf.setRangePosition(vort)
    I(:, i) = vf.impulse(center, 0);
    vfs{i} = vf;
end

ii = sqrt(sum(I.^2, 1));

% Add noise.
u_mean = vfs{1}.meanSpeed(0, 0);
noise_props = 2;
props_count = size(noise_props, 2);

% Number of iterations.
num_ite = 200;

% Impulse values computed with noise.
In = zeros(3, num_frame, props_count, num_ite);

% Unnormalized error of noisy impulse.
err = zeros(3, num_frame, props_count, num_ite);

for k = 1: num_ite
    for p = 1: props_count
        % Compute impulse with noise.
        for i = 1: num_frame
            vf = vfs{i};
            % Introduce noise.
            vf.clearNoise();
            vf.noise_uniform(noise_props(p)*u_mean);
            In(:, i, p, k) = vf.impulse(center, 1);
            err(:, i, p, k) = abs(vf.impulse(center, 1) - I(:, i));
        end
    end
end

% Plot dimensionally.
dims = [2];

for dim = dims
    figure;
    In_dim = squeeze(In(dim,:,:,:));
    err_dim = squeeze(err(dim,:,:,:));
    
    % Average across trials and error levels.
    errorbar(t, mean(In_dim, [2 3]), std(In_dim, 0, [2 3]), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
%     % Average across trials and error levels.
%     errorbar(t, mean(err_dim, [2 3]), std(err_dim, 0, [2 3]), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    
    legend({'noisy'})
    xlabel('$t$')
    ylabel(strcat('$I_', vfs{1}.dim_str(dim), '(t)$'))
    title(strcat('$', vfs{1}.dim_str(dim), '$ Impulse with noise and smoothing'))
end

mag_err = sqrt(sum(err.^2, 1));
mean_mag_err = squeeze(mean(mag_err, 4));
dev_mag_err = squeeze(std(mag_err, 0, 4));

figure;
errorbar(t, mean_mag_err, dev_mag_err, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
xlabel('$t$')
ylabel('$\delta I$')
title('Error Magnitude')