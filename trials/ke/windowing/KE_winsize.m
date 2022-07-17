function [dK, dK_box, dK_gss, dK0, bias_box, bias_gss, dK_sd, dK_sd_box, dK_sd_gss, axes] = ...
    KE_winsize(vf, K0, KEf, props, winsizes, overlap, num_ite, display_plots)
% Vary the window size used in downsampling and present its effect on
% error. The errors returned here are not absolute but the sign is
% retained.

if ~isvector(winsizes)
    error('Uniform windows expected!')
end

% Assume uniform windows in x, y, z.
windows_count = length(winsizes);

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
            KE_err_run_constN(vf, props, K0, KEf, num_ite, [winsizes(k), overlap], {});
end

%%%%%%%%%%% Error plots %%%%%%%%%%%%%%
if ~exist('display_plots', 'var') || ((isinteger(display_plots)||islogical(display_plots)) && ~display_plots)
    return
end

% Handle for figures to return.
axes = {};

% Plot of baseline resolution error.
if ismember('win', display_plots)
    % figure;
    axes{end+1} = scatter(winsizes, dK0, 'k', 'filled');
    xticks(winsizes)
    xlabel('Window size')
    ylabel('$\frac{\delta K}{K}$', 'HorizontalAlignment', 'right', 'Rotation', 0)
    title('Kinetic energy windowing resolution error')
    % Log plot.
    ax = gca;
    ax.XScale = 'log';
    xlim([winsizes(1)/2 winsizes(end)*2])
end

% Plot of windowing error filter errors.
if ismember('resol', display_plots)
    %             figure;
    axes{end+1} = scatter(winsizes, dK0, 'k', 'filled');
%     hold on
%     scatter(winsizes, bias_box, 'r', 'filled', 'Marker', 's')
    hold on
    scatter(winsizes, bias_gss, 'b', 'filled', 'Marker', '^')
    
    xticks(winsizes)
    legend({'unfiltered', 'Gaussian'})
    %             legend({'unfiltered', 'box', 'Gaussian'})
    xlabel('Window size')
    ylabel('$\frac{\delta K}{K}$', 'HorizontalAlignment', 'right', 'Rotation', 0)
    title('Kinetic energy windowing resolution error')
    % Log plot.
    ax = gca;
    ax.XScale = 'log';
    xlim([winsizes(1)/2 winsizes(end)*2])
end

% % Plot of mean errors.
% if plot_signed
%     figure;
%     errorbar(windows, dK, dK_sd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
%     hold on
%     errorbar(windows, dK_box, dK_sd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
%     hold on
%     errorbar(windows, dK_gss, dK_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
%     
%     xticks(windows)
%     legend({'unfiltered', 'box', 'Gaussian'})
%     xlabel('Window size')
%     ylabel('Normalized error')
%     title('KE error with noise')
% end
% 
% if plot_abs
%     figure;
%     errorbar(windows, abs(dK), dK_sd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
%     hold on
%     errorbar(windows, abs(dK_box), dK_sd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
%     hold on
%     errorbar(windows, abs(dK_gss), dK_sd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
%     
%     xticks(windows)
%     legend({'unfiltered', 'box', 'Gaussian'})
%     xlabel('Window size')
%     ylabel('Normalized error')
%     title('Absolute KE error with noise')
% end
