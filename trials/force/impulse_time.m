% Compute the impulse and its derivetive over time without adding noise.

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
% Fixed range where impulse is calculated.
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

% Smooth impulse values temporally and then compute time-derivative.
[I_fit, dI_dt, spls] = smoothVector_temporal(I, t);

% Magnitude of impulse.
i = sqrt(sum(I_fit.^2, 1));
% This is not the derivative of magnitude, but the magnitude of derivative.
di_dt = sqrt(sum(dI_dt.^2, 1));

% Plot impulse and dI/dt.
% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [2];
for dim = dims
    % Without smoothing.
    Idim = I(dim, :);
    
    figure;
    scatter(t, Idim, 'filled', 'b')
    hold on
    fnplt(spls{dim})
    
    legend({'unfiltered', 'spline smoothed'})
    xlabel('$t$')
    ylabel(strcat('$I_', vfs{1}.dim_str(dim), '$'))
    title('Fixed Region and Stationary Origin')
    
    figure;
    scatter(t, dI_dt(dim,:), 'filled', 'r')
    xlabel('$t$')
    ylabel(strcat('$\frac{\partial I_', vfs{1}.dim_str(dim), '}{\partial t}$'))
    title('Fixed Region and Stationary Origin')
end

% Hereafter the spline-fitted impulse values and first derivative computed
% above are taken as correct.
u_mean = vfs{1}.meanSpeed(0, 0);
noise_props = 0: 0.5: 2;
props_count = size(noise_props, 2);

% Impulse values computed with noise.
In = zeros(3, num_frame, props_count);
% Noisy impulse value with smoothing.
I_fitn = zeros(3, num_frame, props_count);
% Derivatives computed with noise and smoothing.
dIn_dt = zeros(3, num_frame, props_count);
% Deviation of impulse value given noise.
err_I = zeros(props_count, 2);
% Deviation of temporal derivative given noise.
err_dI = zeros(props_count, 1);

% Cell array for smoothing splines fitted.
spls = cell(3, props_count);

for p = 1: props_count
    % Compute impulse with noise.
    for i = 1: num_frame
        vf = vfs{i};
        % Introduce noise.
        vf.clearNoise();
        vf.noise_uniform(noise_props(p)*u_mean);
        In(:, i, p) = vf.impulse(center, 1);
    end
    
    % Mean percentage of error without smoothing.
    err_I(p, 1) = mean(sqrt(sum((In(:,:,p) - I_fit).^2, 1)) ./ i);
    
    % Apply smoothing splines dimensionally.
    [I_s, dI_s, sps] = smoothVector_temporal(In(:,:,p), t);
    I_fitn(:,:,p) = I_s;
    dIn_dt(:,:,p) = dI_s;
    % Indexing convention for cell arrays.
    spls{1, p} = sps{1};
    spls{2, p} = sps{2};
    spls{3, p} = sps{3};
    
    % Error after smoothing.
    err_I(p, 2) = mean(sqrt(sum((I_s - I_fit).^2, 1)) ./ i);
    
    % Error of derivative.
    err_dI(p) = mean(sqrt(sum((dI_s - dI_dt).^2, 1)) ./ di_dt);
end

% Plot dimensionally.
dims = [2];

for dim = dims
    figure;
    In_dim = squeeze(In(dim,:,:));
    % Noisy impulse values.
    errorbar(t, mean(In_dim, 2), std(In_dim, 1, 2), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    
    % Plot fitted splines.
    for p = 1: props_count
        hold on
        fnplt(spls{dim, p})
    end
    
    legend([{'noisy'}; split(cellstr(num2str(noise_props)))])
    xlabel('$t$')
    ylabel(strcat('$I_', vfs{1}.dim_str(dim), '(t)$'))
    title(strcat('$', vfs{1}.dim_str(dim), '$ Impulse with noise and smoothing'))
end

% Plot error with respect to noise.
figure;
scatter(noise_props, err_I(:, 1), 'filled')
hold on
scatter(noise_props, err_I(:, 2), 'filled')

legend('noisy', 'smoothed')
xlabel('$\delta u$')
ylabel('$|\delta I|$')
title('Normalized impulse error averaged over time')

