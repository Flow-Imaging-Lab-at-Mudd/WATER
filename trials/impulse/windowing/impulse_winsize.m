function [dI, dI_box, dI_gss, dI0, bias_box, bias_gss, ...
    di, di_box, di_gss, di0, mag_bias_box, mag_bias_gss, ...
    dI_sd, dI_sd_box, dI_sd_gss, di_sd, di_sd_box, di_sd_gss] = ...
    impulse_winsize(vf, I0, origin, props, winsizes, overlap, display_plots)
% Vary the window size used in downsampling and present its effect on error.
%
% April, 2022

if ~isvector(winsizes)
    error('Uniform windows expected!')
end
% Assume uniform windows in x, y, z.
win_count = length(winsizes);

% Number of repetitions for noise trial.
num_ite = 10;

% Containers for error data at different window sizes.
dI = zeros(3, win_count);
dI_sd = zeros(3, win_count);
dI0 = zeros(3, win_count);
dI_box = zeros(3, win_count);
dI_sd_box = zeros(3, win_count);
dI_gss = zeros(3, win_count);
dI_sd_gss = zeros(3, win_count);
bias_box = zeros(3, win_count);
bias_gss = zeros(3, win_count);

% Magnitude errors.
di = zeros(1, win_count);
di_sd = zeros(1, win_count);
di0 = zeros(1, win_count);
di_box = zeros(1, win_count);
di_sd_box = zeros(1, win_count);
mag_bias_box = zeros(1, win_count);
di_gss = zeros(1, win_count);
di_sd_gss = zeros(1, win_count);
mag_bias_gss = zeros(1, win_count);

for k = 1: win_count
    [dI(:,k), dI_box(:,k), dI_gss(:,k), dI0(:,k), bias_box(:,k), bias_gss(:,k), ...
        dI_sd(:,k), dI_sd_box(:,k), dI_sd_gss(:,k), ...
        di(:,k), di_box(:,k), di_gss(:,k), di0(:,k), mag_bias_box(:,k), mag_bias_gss(:,k), ...
        di_sd(:,k), di_sd_box(:,k), di_sd_gss(:,k), ~] = ...
        impulse_err_run_constN(vf, props, origin, I0, num_ite, [winsizes(k), overlap]);
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
    scatter(winsizes, dI0(dim,:), 'k', 'filled')
    xticks(winsizes)
    xlabel('Window size')
    ylabel('Normalized error')
    title(sprintf('Windowing resolution error in $\\hat{%s}$', dim_str{dim}))
    
    % Plot of filter errors.
    figure;
    scatter(winsizes, dI0(dim,:), 'k', 'filled')
    hold on
    scatter(winsizes, bias_box(dim,:), 'r', 'filled')
    hold on
    scatter(winsizes, bias_gss(dim,:), 'b', 'filled')
    
    xticks(winsizes)
    legend({'windowing', 'box', 'Gaussian'})
    xlabel('Window size')
    ylabel('Normalized error')
    title(sprintf('Filter errors in $\\hat{%s}$', dim_str{dim}))
    
    % Plot of mean errors.
    figure;
    errorbar(winsizes, dI(dim,:), dI_sd(dim,:), 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    errorbar(winsizes, dI_box(dim,:), dI_sd_box(dim,:), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    errorbar(winsizes, dI_gss(dim,:), dI_sd_gss(dim,:), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
    
    xticks(winsizes)
    legend({'unfiltered', 'box', 'Gaussian'})
    xlabel('Window size')
    ylabel('Normalized error')
    title(sprintf('Error in $\\hat{%s}$ with noise', dim_str{dim}))
end

%%%%%%%%%%% Magnitude plots %%%%%%%%%%%%%%
% Plot of windowing error.
figure;
scatter(winsizes, di0, 'k', 'filled')
xticks(winsizes)
xlabel('Window size')
ylabel('Normalized error')
title('Windowing resolution error of impulse')

% Plot of bias.
figure;
scatter(winsizes, di0, 'k', 'filled')
hold on
scatter(winsizes, mag_bias_box, 'r', 'filled')
hold on
scatter(winsizes, mag_bias_gss, 'b', 'filled')

xticks(winsizes)
legend({'windowing', 'box', 'Gaussian'})
xlabel('Window size')
ylabel('Normalized error')
title('Baseline impulse resolution error of filters')

% Plot of mean absolute errors.
figure;
errorbar(winsizes, di, di_sd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
hold on
errorbar(winsizes, di_box, di_sd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
hold on
errorbar(winsizes, di_gss, di_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)

xticks(winsizes)
legend({'unfiltered', 'box', 'Gaussian'})
xlabel('Window size')
ylabel('Normalized error')
title('Impulse error with noise')
