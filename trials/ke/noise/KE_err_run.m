function [dK, dK_box, dK_gss, bias_box, bias_gss] = KE_err_run(vf, props, K0, KEf, ...
    display_plots)
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
% 'display_plots' is a boolean for displaying plots generated noise
% propagation plots herein. Which specific types of plots are chosen must
% be indicated inside this file. Magnitude plots are by deafult displayed
% if a true value is passed in.
%
% Derek Li, March, 2022

% Set constant maximal magnitude of noise.
u_mean = vf.meanSpeed(0, 0);

% Error in energy estimation given noise.
dK = zeros(size(props));
% Box smoothing.
dK_box = zeros(size(props));
% Gaussian smoothing.
dK_gss = zeros(size(props));

% Plot energy estimation error for small and large values of noise.
for i = 1: size(props, 2)
    vf.clearNoise();
    N = vf.noise_uniform(props(i)*u_mean);
    dK(i) = KEf(vf, 1) - K0;
    % Result with box smoothing.
    vf.smoothNoise('box');
    dK_box(i) = KEf(vf, 1) - K0;
    % Reset and smooth with gaussian filter.
    vf.setNoise(N)
    vf.smoothNoise('gaussian');
    dK_gss(i) = KEf(vf, 1) - K0;
end

% Normalize.
dK = dK / K0;
dK_box = dK_box / K0;
dK_gss = dK_gss / K0;

% Baseline smoother biases.
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
    scatter(props, dK, 'k', 'filled')
    hold on
    scatter(props, dK_box, 'r', 'filled')
    hold on
    scatter(props, dK_gss, 'b', 'filled')
    legend({'unfiltered', 'box-filtered', 'Gaussian-filtered'})
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    ylabel('$\frac{\delta K}{K}$')
    title('Normalized KE Error')
    ax = gca;
    ax.FontSize = fsize;
end

%%%%%%%%%%%%%%%%%% Plot absolute KE error %%%%%%%%%%%%%%%%%%%%%
plot_abs_err = 0;
if plot_abs_err
    abs_dK = abs(dK);
    abs_dK_box = abs(dK_box);
    abs_dK_gss = abs(dK_gss);
    figure;
    scatter(props, abs_dK, 'k', 'filled')
    hold on
    scatter(props, abs_dK_box, 'r', 'filled')
    hold on
    scatter(props, abs_dK_gss, 'b', 'filled')
    
    legend({'unfiltered', 'box-filtered', 'Gaussian-filtered'}, 'Interpreter', 'latex')
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    ylabel('$|\frac{\delta K}{K}|$')
    title('Absolute KE Error')
    ax = gca;
    ax.FontSize = fsize;
end