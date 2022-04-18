function [dK, dK_box, dK_gss, dK0, bias_box, bias_gss, dK_sd, dK_sd_box, dK_sd_gss] = ...
    KE_overlap(vf, K0, KEf, props, window, overlaps, display_plots)
% Vary the overlap used in downsampling and present its effect on error.

if ~isvector(overlaps)
    error('Uniform overlaps ratios expected!')
end

% Assume uniform windows in x, y, z.
ops_count = length(overlaps);

% Number of repetitions for noise trial.
num_ite = 10;

% Errors with noise.
dK = zeros(ops_count, 1);
dK_sd = zeros(ops_count, 1);
dK_box = zeros(ops_count, 1);
dK_sd_box = zeros(ops_count, 1);
dK_gss = zeros(ops_count, 1);
dK_sd_gss = zeros(ops_count, 1);

% Resolution errors.
bias_box = zeros(ops_count, 1);
bias_gss = zeros(ops_count, 1);
dK0 = zeros(ops_count, 1);

for k = 1: ops_count
    [dK(k), dK_box(k), dK_gss(k), dK0(k), bias_box(k), bias_gss(k), ...
        dK_sd(k), dK_sd_box(k), dK_sd_gss(k), vfd] = ...
            KE_err_run_constN(vf, props, K0, KEf, num_ite, [window, overlaps(k)]);
        disp(vfd.dims)
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
    scatter(overlaps, dK0, 'k', 'filled')
    xticks(overlaps)
    xlabel('Overlap ratio of windowing')
    ylabel('Normalized error')
    title('Windowing resolution error of KE')
end

if plot_abs
    figure;
    scatter(overlaps, abs(dK0), 'k', 'filled')
    xticks(overlaps)
    xlabel('Overlap ratio of windowing')
    ylabel('Normalized error')
    title('Absolute windowing resolution error of KE')
end

% Plot of bias.
if plot_signed
    figure;
    scatter(overlaps, dK0, 'k', 'filled')
    hold on
    scatter(overlaps, bias_box, 'r', 'filled')
    hold on
    scatter(overlaps, bias_gss, 'b', 'filled')
    xticks(overlaps)
    legend({'windowing', 'box', 'Gaussian'})
    xlabel('Overlap ratio of windowing')
    ylabel('Normalized error')
    title('Baseline KE resolution error of filters')
end

if plot_abs
    figure;
    scatter(overlaps, abs(dK0), 'k', 'filled')
    hold on
    scatter(overlaps, abs(bias_box), 'r', 'filled')
    hold on
    scatter(overlaps, abs(bias_gss), 'b', 'filled')
    xticks(overlaps)
    legend({'windowing', 'box', 'Gaussian'})
    xlabel('Overlap ratio of windowing')
    ylabel('Normalized error')
    title('Absolute baseline KE resolution error of filters')
end

% Plot of mean absolute errors.
if plot_signed
    figure;
    errorbar(overlaps, dK, dK_sd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    errorbar(overlaps, dK_box, dK_sd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    errorbar(overlaps, dK_gss, dK_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
    xticks(overlaps)
    legend({'unfiltered', 'box', 'Gaussian'})
    xlabel('Overlap ratio of windowing')
    ylabel('Normalized error')
    title('KE Error with noise')
end

if plot_abs
    figure;
    errorbar(overlaps, abs(dK), dK_sd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    errorbar(overlaps, abs(dK_box), dK_sd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    errorbar(overlaps, abs(dK_gss), dK_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1) 
    xticks(overlaps)
    legend({'unfiltered', 'box', 'Gaussian'})
    xlabel('Overlap ratio of windowing')
    ylabel('Normalized error')
    title('KE Error with noise')
end

