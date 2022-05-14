function [di, di_box, di_gss, di0, mag_bias_box, mag_bias_gss, di_sd, di_box_sd, di_gss_sd] = ...
    impulse_winsize(vf, I0, origin, props, windows, overlap, display_plots)
% Vary the window size used in downsampling and present its effect on error.
%
% April, 2022

if ~isvector(windows)
    error('Uniform windows expected!')
end
% Assume uniform windows in x, y, z.
windows_count = length(windows);

% Number of repetitions for noise trial.
num_ite = 5;

% Containers for error data at different window sizes.
dI = zeros(3, windows_count, num_ite);
dI0 = zeros(3, windows_count);
dI_box = zeros(3, windows_count, num_ite);
dI_gss = zeros(3, windows_count, num_ite);
bias_box = zeros(3, windows_count);
bias_gss = zeros(3, windows_count);

for k = 1: windows_count
    vfd = vf.downsample(windows(k), overlap, 0);
    for i = 1: num_ite
        % Run script for impulse error sampling.
        [dI(:,k,i), ~, dI_box(:,k,i), ~, dI_gss(:,k,i), ~, bias_box(:,k), ...
            bias_gss(:,k), dI0(:,k)] = impulse_err_stats(vfd, props, origin, I0, []);
    end
end

% Magnitudes.
di = squeeze(sqrt(sum(dI.^2, 1)));
di_box = squeeze(sqrt(sum(dI_box.^2, 1)));
di_gss = squeeze(sqrt(sum(dI_gss.^2, 1)));

% Average over trials.
di_sd = std(di, 0, 2);
di_box_sd = std(di_box, 0, 2);
di_gss_sd = std(di_gss, 0, 2);
dI_sd = std(dI, 0, 3);
dI_box_sd = std(dI_box, 0, 3);
dI_gss_sd = std(dI_gss, 0, 3);

dI = mean(dI, 3);
dI_box = mean(dI_box, 3);
dI_gss = mean(dI_gss, 3);
di = mean(di, 2);
di_box = mean(di_box, 2);
di_gss = mean(di_gss, 2);

di0 = sqrt(sum(dI0.^2, 1));
mag_bias_box = sqrt(sum(bias_box.^2, 1));
mag_bias_gss = sqrt(sum(bias_gss.^2, 1));

%%%%%%%%%%% Dimensional Plots of Error %%%%%%%%%%%%%
if ~exist('display_plots', 'var') || ~display_plots
    return
end
% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [];
dim_str = {'x', 'y', 'z'};

for dim = dims
    % Plot of baseline resolution error.
    figure;
    scatter(windows, dI0(dim,:), 'k', 'filled')
    xticks(windows)
    xlabel('Window size')
    ylabel('Normalized error')
    title(sprintf('Windowing resolution error in $\\hat{%s}$', dim_str{dim}))
    
    % Plot of filter errors.
    figure;
    scatter(windows, dI0(dim,:), 'k', 'filled')
    hold on
    scatter(windows, bias_box(dim,:), 'r', 'filled')
    hold on
    scatter(windows, bias_gss(dim,:), 'b', 'filled')
    
    xticks(windows)
    legend({'windowing', 'box', 'Gaussian'})
    xlabel('Window size')
    ylabel('Normalized error')
    title(sprintf('Filter errors in $\\hat{%s}$', dim_str{dim}))
    
    % Plot of mean errors.
    figure;
    errorbar(windows, dI(dim,:), dI_sd(dim,:), 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    errorbar(windows, dI_box(dim,:), dI_box_sd(dim,:), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    errorbar(windows, dI_gss(dim,:), dI_gss_sd(dim,:), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
    
    xticks(windows)
    legend({'unfiltered', 'box', 'Gaussian'})
    xlabel('Window size')
    ylabel('Normalized error')
    title(sprintf('Error in $\\hat{%s}$ with noise', dim_str{dim}))
end

%%%%%%%%%%% Magnitude plots %%%%%%%%%%%%%%
% Plot of windowing error.
figure;
scatter(windows, di0, 'k', 'filled')
xticks(windows)
xlabel('Window size')
ylabel('Normalized error')
title('Windowing resolution error of impulse')

% Plot of bias.
figure;
scatter(windows, di0, 'k', 'filled')
hold on
scatter(windows, mag_bias_box, 'r', 'filled')
hold on
scatter(windows, mag_bias_gss, 'b', 'filled')

xticks(windows)
legend({'windowing', 'box', 'Gaussian'})
xlabel('Window size')
ylabel('Normalized error')
title('Baseline impulse resolution error of filters')

% Plot of mean absolute errors.
figure;
errorbar(windows, di, di_sd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
hold on
errorbar(windows, di_box, di_box_sd, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
errorbar(windows, di_gss, di_gss_sd, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)

xticks(windows)
legend({'unfiltered', 'box', 'Gaussian'})
xlabel('Window size')
ylabel('Normalized error')
title('Impulse error with noise')
