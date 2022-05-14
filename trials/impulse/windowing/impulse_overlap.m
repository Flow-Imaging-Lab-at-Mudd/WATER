function [dI, dI_box, dI_gss, dI0, bias_box, bias_gss, ...
    di, di_box, di_gss, di0, mag_bias_box, mag_bias_gss, ...
    dI_sd, dI_sd_box, dI_sd_gss, di_sd, di_sd_box, di_sd_gss] = ...
    impulse_overlap(vf, I0, origin, props, winsize, overlaps, display_plots)
% Vary the overlap ratio used in downsampling and present its effect on error.
%
% April, 2022

if ~isvector(overlaps)
    error('Uniform overlaps ratios expected!')
end

% Assume uniform overlap ratios in x, y, z.
ops_count = length(overlaps);

% Number of repetitions for noise trial.
num_ite = 10;

% Containers for error data at different window sizes.
dI = zeros(3, ops_count);
dI_sd = zeros(3, ops_count);
dI0 = zeros(3, ops_count);
dI_box = zeros(3, ops_count);
dI_sd_box = zeros(3, ops_count);
dI_gss = zeros(3, ops_count);
dI_sd_gss = zeros(3, ops_count);
bias_box = zeros(3, ops_count);
bias_gss = zeros(3, ops_count);

% Magnitude errors.
di = zeros(1, ops_count);
di_sd = zeros(1, ops_count);
di0 = zeros(1, ops_count);
di_box = zeros(1, ops_count);
di_sd_box = zeros(1, ops_count);
mag_bias_box = zeros(1, ops_count);
di_gss = zeros(1, ops_count);
di_sd_gss = zeros(1, ops_count);
mag_bias_gss = zeros(1, ops_count);

for k = 1: ops_count
    [dI(:,k), dI_box(:,k), dI_gss(:,k), dI0(:,k), bias_box(:,k), bias_gss(:,k), ...
        dI_sd(:,k), dI_sd_box(:,k), dI_sd_gss(:,k), ...
        di(:,k), di_box(:,k), di_gss(:,k), di0(:,k), mag_bias_box(:,k), mag_bias_gss(:,k), ...
        di_sd(:,k), di_sd_box(:,k), di_sd_gss(:,k), ~] = ...
        impulse_err_run_constN(vf, props, origin, I0, num_ite, [winsize, overlaps(k)]);
end

%%%%%%%%%%% Dimensional Plots of Error %%%%%%%%%%%%%
if ~exist('display_plots', 'var') || ~display_plots
    return
end
% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [3];
dim_str = {'x', 'y', 'z'};

for dim = dims
    % Plot of baseline resolution error.
    figure;
    scatter(overlaps, dI0(dim,:), 'k', 'filled')
    xticks(overlaps)
    xlabel('Overlap ratio of window')
    ylabel('Normalized error')
    title(sprintf('Windowing resolution error in $\\hat{%s}$', dim_str{dim}))
    
    % Plot of filter errors.
    figure;
    scatter(overlaps, dI0(dim,:), 'k', 'filled')
    hold on
    scatter(overlaps, bias_box(dim,:), 'r', 'filled')
    hold on
    scatter(overlaps, bias_gss(dim,:), 'b', 'filled')
    
    xticks(overlaps)
    legend({'windowing', 'box', 'Gaussian'})
    xlabel('Overlap ratio of window')
    ylabel('Normalized error')
    title(sprintf('Filter errors in $\\hat{%s}$', dim_str{dim}))
    
    % Plot of mean errors.
    figure;
    errorbar(overlaps, dI(dim,:), dI_sd(dim,:), 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    errorbar(overlaps, dI_box(dim,:), dI_sd_box(dim,:), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    errorbar(overlaps, dI_gss(dim,:), dI_sd_gss(dim,:), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
    
    xticks(overlaps)
    legend({'unfiltered', 'box', 'Gaussian'})
    xlabel('Overlap ratio of window')
    ylabel('Normalized error')
    title(sprintf('Error in $\\hat{%s}$ with noise', dim_str{dim}))
end

%%%%%%%%%%% Magnitude plots %%%%%%%%%%%%%%
% Plot of windowing error.
figure;
scatter(overlaps, di0, 'k', 'filled')
xticks(overlaps)
xlabel('Overlap ratio of window')
ylabel('Normalized error')
title('Windowing resolution error of impulse')

% Plot of bias.
figure;
scatter(overlaps, di0, 'k', 'filled')
hold on
scatter(overlaps, mag_bias_box, 'r', 'filled')
hold on
scatter(overlaps, mag_bias_gss, 'b', 'filled')

xticks(overlaps)
legend({'windowing', 'box', 'Gaussian'})
xlabel('Overlap ratio of window')
ylabel('Normalized error')
title('Baseline impulse resolution error of filters')

% Plot of mean absolute errors.
figure;
errorbar(overlaps, di, di_sd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
hold on
errorbar(overlaps, di_box, di_sd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
errorbar(overlaps, di_gss, di_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)

xticks(overlaps)
legend({'unfiltered', 'box', 'Gaussian'})
xlabel('Overlap ratio of window')
ylabel('Normalized error')
title('Impulse error with noise')
