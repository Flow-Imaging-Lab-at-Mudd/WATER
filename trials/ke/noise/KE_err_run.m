function [dK, dK_box, dK_gss, dK0, bias_box, bias_gss, dKsd, dKsd_box, dKsd_gss, vf] = ...
    KE_err_run(vf, props, K0, KEf, num_ite, window_params, display_plots)
% An noise-propagation trial of computing KE after introducing specified
% levels of noise.
% 
% The theoretical (expected) kinetic energy of the currently effective
% region is provided as 'K0' to determine the error.
% 
% A customized function handle for computing kinetic energy is given by
% 'KEf', which is of form KEf(vf, with_noise), where the first input is a
% VF object and the second a boolean for whether noise is to be included.
% Such a handle is useful, for example, in selecting only a portion, say
% the vortical, of the velocity field to compute KE.
%
% Introduce levels of noise proportional to the mean speed in the effective
% region, according to 'props', e.g. 0: 0.1: 3. It is required that
% props(1) = 0 so that baseline resolution errors can be computed.
%
% 'vf' is presume to have range properly set.
%
% 'num_ite' specifies the number of iterations the computation is to be
% repeated and their results averaged to account for stochasticity.
%
% 'display_plots' is a boolean for displaying plots generated noise
% propagation plots herein. Which specific types of plots are chosen must
% be indicated inside this file. Magnitude plots are by deafult displayed
% if a true value is passed in.
%
% Derek Li, March, 2022

% Optional windowing operation.
if isvector(window_params) && length(window_params) == 2
    winsize = window_params(1);
    overlap = window_params(2);
    vf = vf.downsample(winsize, overlap, 0);
end

% Set constant maximal magnitude of noise.
u_mean = vf.meanSpeed(0, 0);

% Error in energy estimation given noise.
dK = zeros(length(props), num_ite);
% Box smoothing.
dK_box = zeros(length(props), num_ite);
% Gaussian smoothing.
dK_gss = zeros(length(props), num_ite);

% Plot energy estimation error for small and large values of noise.
for i = 1: size(props, 2)
    for j = 1: num_ite
        vf.clearNoise();
        N = vf.noise_uniform(props(i)*u_mean);
        dK(i,j) = KEf(vf, 1) - K0;
        % Result with box smoothing.
        vf.smoothNoise('box');
        dK_box(i,j) = KEf(vf, 1) - K0;
        % Reset and smooth with gaussian filter.
        vf.setNoise(N)
        vf.smoothNoise('gaussian');
        dK_gss(i,j) = KEf(vf, 1) - K0;
    end
end

% Normalize.
dK = dK / K0;
dK_box = dK_box / K0;
dK_gss = dK_gss / K0;

% Average across trials.
dKsd = std(dK, 0, 2);
dKsd_box = std(dK_box, 0, 2);
dKsd_gss = std(dK_gss, 0, 2);

dK = mean(dK, 2);
dK_box = mean(dK_box, 2);
dK_gss = mean(dK_gss, 2);

% Baseline resolution errors.
dK0 = dK(1);
bias_box = dK_box(1);
bias_gss = dK_gss(1);

%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%
% Font size for plots.
fsize = 10;

if ~exist('display_plots', 'var') || ~display_plots
    return
end

plot_signed_err = 1;
if plot_signed_err
    figure;
    errorbar(props, dK, dKsd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    errorbar(props, dK_box, dKsd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    errorbar(props, dK_gss, dKsd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
    
    legend({'unfiltered', 'box-filtered', 'Gaussian-filtered'})
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    ylabel('$\frac{\delta K}{K}$')
    title('Normalized KE Error')
    ax = gca;
    ax.FontSize = fsize;
end

%%%%%%%%%%%%%%%%%% Plot absolute KE error %%%%%%%%%%%%%%%%%%%%%
plot_abs_err = 1;
if plot_abs_err
    abs_dK = abs(dK);
    abs_dK_box = abs(dK_box);
    abs_dK_gss = abs(dK_gss);
    abs_dKd = abs(dKsd);
    abs_dKd_box = abs(dKsd_box);
    abs_dKd_gss = abs(dKsd_gss);
    
    figure;
    errorbar(props, abs_dK, abs_dKd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
    hold on
    errorbar(props, abs_dK_box, abs_dKd_box, 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
    hold on
    errorbar(props, abs_dK_gss, abs_dKd_gss, 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
    
    legend({'unfiltered', 'box-filtered', 'Gaussian-filtered'}, 'Interpreter', 'latex')
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    ylabel('$|\frac{\delta K}{K}|$')
    title('Absolute KE Error')
    ax = gca;
    ax.FontSize = fsize;
end