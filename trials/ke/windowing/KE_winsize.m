function [dK, dK_box, dK_gss, dK0, bias_box, bias_gss, dK_sd, dK_sd_box, dK_sd_gss] = ...
    KE_winsize(vf, K0, KEf, props, windows, overlap, display_plots)
% Vary the window size used in downsampling and present its effect on
% error. The errors returned here are not absolute but the sign is
% retained.

if ~isvector(windows)
    error('Uniform windows expected!')
end

% Assume uniform windows in x, y, z.
windows_count = length(windows);

% Number of repetitions for noise trial.
num_ite = 10;

% Errors with noise.
dK = zeros(windows_count, 1);
dK_sd = zeros(windows_count, 1);
dK_box = zeros(windows_count, 1);
dK_sd_box = zeros(windows_count, 1);
dK_gss = zeros(windows_count, 1);
dK_sd_gss = zeros(windows_count, 1);

% Resolution errors.
bias_box = zeros(windows_count, 1);
bias_gss = zeros(windows_count, 1);
dK0 = zeros(windows_count, 1);

for k = 1: windows_count
    [dK(k), dK_box(k), dK_gss(k), dK0(k), bias_box(k), bias_gss(k), ...
        dK_sd(k), dK_sd_box(k), dK_sd_gss(k)] = ...
            KE_err_run_constN(vf, props, K0, KEf, num_ite, [windows(k), overlap]);
end

%%%%%%%%%%% Error plots %%%%%%%%%%%%%%
if ~exist('display_plots', 'var') || ~display_plots
    return
end

% Whether error is to be plotted with sign (can be both).
plot_signed = 1;
plot_abs = 0;

 % Plot of windowing error.
if plot_signed
    figure;
    scatter(windows, dK0, 'k', 'filled')
    xticks(windows)
    xlabel('Window size')
    ylabel('Normalized error')
    title('Windowing resolution error of KE')
end

if plot_abs
    figure;
    scatter(windows, abs(dK0), 'k', 'filled')
    xticks(windows)
    xlabel('Window size')
    ylabel('Normalized error')
    title('Absolute windowing resolution error of KE')
end

% Plot of bias.
if plot_signed
    figure;
    scatter(windows, dK0, 'k', 'filled')
    hold on
    scatter(windows, bias_box, 'r', 'filled')
    hold on
    scatter(windows, bias_gss, 'b', 'filled')
    
    xticks(windows)
    legend({'windowing', 'box', 'Gaussian'})
    xlabel('Window size')
    ylabel('Normalized error')
    title('Baseline KE resolution error of filters')
end

if plot_abs
    figure;
    scatter(windows, abs(dK0), 'k', 'filled')
    hold on
    scatter(windows, abs(bias_box), 'r', 'filled')
    hold on
    scatter(windows, abs(bias_gss), 'b', 'filled')
    
    xticks(windows)
    legend({'windowing', 'box', 'Gaussian'})
    xlabel('Window size')
    ylabel('Normalized error')
    title('Absolute baseline KE resolution error of filters')
end

% Plot of mean errors.
if plot_signed
    figure;
    errorbar(windows, dK, dK_sd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    errorbar(windows, dK_box, dK_sd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    errorbar(windows, dK_gss, dK_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
    
    xticks(windows)
    legend({'unfiltered', 'box', 'Gaussian'})
    xlabel('Window size')
    ylabel('Normalized error')
    title('KE error with noise')
end

if plot_abs
    figure;
    errorbar(windows, abs(dK), dK_sd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    errorbar(windows, abs(dK_box), dK_sd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    errorbar(windows, abs(dK_gss), dK_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
    
    xticks(windows)
    legend({'unfiltered', 'box', 'Gaussian'})
    xlabel('Window size')
    ylabel('Normalized error')
    title('Absolute KE error with noise')
end
