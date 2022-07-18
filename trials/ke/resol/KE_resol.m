function [fres, dK, dK_box, dK_gss, dK0, bias_box, bias_gss, dK_sd, dK_sd_box, dK_sd_gss, vfds] = ...
    KE_resol(l, vr, u0, min_fres, max_fres, fres_inc, props, err_level, ...
        num_ite, window_params, display_plots, vfds)
% Compute the signed errors of KE computations at various resolutions as
% specified with smoothing options.

% Radius of vortex.
fr = l*vr;
% Desired feature resolutions.
fres = min_fres: fres_inc: max_fres;

% If a cell array of VelocityField objects are given, assume they are of
% the appropriate resolutions and already downsampled.
if ~exist('vfds', 'var')
    vfds = cell(1, length(fres));
else
    % Non-trivial windowing parameters will be ignored if an array of VFs
    % is provided.
    window_params = [];
end

% Configure windowing parameters, if applicable.
windowing = false;
if isvector(window_params) && length(window_params) == 2
    windowing = true;
    winsize = window_params(1);
    overlap = window_params(2);
end

% Generate range of spacings for evenly spaced feature resolutions.
if ~windowing
    % Global spacing of downsampled data.
    sps = zeros(1, floor((max_fres-min_fres)/fres_inc) + 1);
    % Minimal feature resolution.
    sps(1) = fr / min_fres;
    for i = 2: size(sps, 2)
        sps(i) = fr / (fres_inc + fr/sps(i-1));
    end
else
    % Compute initial resolutions to produce the desired given feature
    % resolutions after downsampling.
    op = round(winsize .* overlap);
    op = min([op; winsize-1], [], 1);
    sps = 2*fr ./ (2*fres * (winsize-op) + winsize - 1);
end
sps_count = size(sps, 2);
% Normalize spacing.
spr = sps / fr;

% Record feature resolutions.
fres = zeros(1, sps_count);
% Containers for data across all runs.
% Errors here are mean absolute errors.
dK = zeros(sps_count, 1);
dK_sd = zeros(sps_count, 1);
dK_box = zeros(sps_count, 1);
dK_sd_box = zeros(sps_count, 1);
dK_gss = zeros(sps_count, 1);
dK_sd_gss = zeros(sps_count, 1);
% Error magnitudes.
dKm = zeros(sps_count, 1);
dKm_sd = zeros(sps_count, 1);
dKm_box = zeros(sps_count, 1);
dKm_sd_box = zeros(sps_count, 1);
dKm_gss = zeros(sps_count, 1);
dKm_sd_gss = zeros(sps_count, 1);

% Resolution errors.
dK0 = zeros(sps_count, 1);
bias_box = zeros(sps_count, 1);
bias_gss = zeros(sps_count, 1);

% Compute vortical KE of Hill's vortex.
KEf = @(vf, n) Hill_VortKE(vf, fr, n);
% Vortex parameters.
density = 1000;
len_unit = 1e-3;
% Theoretical KE.
K0 = Hill_KE(density, len_unit, fr, u0);

for k = 1: sps_count
    disp(['Spacing ' num2str(k)])
    % Generate VF of appropriate resolution if necessary.
    if isempty(vfds{k})
        % Construct Hill vortex with specified resolution.
        [x, y, z, u, v, w] = Hill_Vortex(spr(k), l, vr, u0, 1);
        vf = VelocityField.importCmps(x, y, z, u, v, w);
        % Focus on vortical region.
        vf.setRangePosition(fr*repmat([-1 1], 3, 1))
    else
        vf = vfds{k};
    end
    
    % Windowing is performed within error trial.
    [dK(k), dK_box(k), dK_gss(k), dK0(k), bias_box(k), bias_gss(k), ...
        dK_sd(k), dK_sd_box(k), dK_sd_gss(k), vfd, ...
        dKm(k), dKm_box(k), dKm_gss(k), dKm_sd(k), dKm_sd_box(k), dKm_sd_gss(k)] = ...
            KE_err_run_constN(vf, props, K0, KEf, num_ite, window_params, {});
    fres(k) = (vfd.span(1)-1)/2;
    % Store this vf to return.
    vfds{k} = vfd;
end

% Compute minimum resolution needed when smoothers are applied. First row
% corresponds to biases; second row, noisy error. Box, first column;
% Gaussian, second.
min_res = -1*ones(2, 2);

% Lowest feature resolutions to achieve the desired error level.
try
    min_res(1,1) = fres(find(abs(dK0) < err_level, 1));
catch
end
try
    min_res(1,2) = fres(find(abs(bias_gss) < err_level, 1));
catch
end

% try
%     min_res(2,1) = fres(find(abs(dKm_box) < err_level, 1));
% catch
% end
try
    min_res(2,2) = fres(find(abs(dKm_gss) < err_level, 1));
catch
end

disp('---------KE minimal resolutions---------')
% Print resolution requirements.
fprintf('Minimum resolution for %.0f%% bias: \n', err_level*100)
fprintf('Original: %d \n', min_res(1,1))
fprintf('Gaussian: %d \n', min_res(1,2))

fprintf('Minimum resolution for %.0f%% noise-propagated error: \n', err_level*100)
% fprintf('Box: %d \n', min_res(2,1))
fprintf('Gaussian: %d \n', min_res(2,2))

%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%
if isempty(display_plots)
    return
end

% Font.
font = 'Arial';
fontSize = 8;

if ismember('bias', display_plots)
    % Smoother bias plot.
%     figure;
    h11=scatter(fres, dK0, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1);
    hold on
    h12=scatter(fres, bias_box, 'r*', 'LineWidth', 1);
    hold on
    h13=scatter(fres, bias_gss, 'b^', 'filled', 'LineWidth', 1);
    xlim([2 28])
    ylim([-0.35 0.05]);
    legend([h11,h12,h13],{'Unfiltered', 'Box', 'Gaussian'},'interpreter','none','location','southeast')
    xlabel('\kappa', 'FontName', font, 'FontSize', 1.25*fontSize,'interpreter','tex')
    ylabel('$\frac{\delta (KE)}{KE}$','FontName',font,'FontSize',1.5*fontSize)
    title('(a) Kinetic energy error','FontName',font,'FontSize',fontSize,'interpreter','latex','fontweight','normal')
    box on
end

if ismember('signed', display_plots)
    % Mean error plot of signed error.
%     figure;
%     errorbar(fres, dK, dK_sd, 'ko', 'MarkerFaceColor','black', 'LineWidth', 1)
%     hold on
    h1=errorbar(fres, dK_box, dK_sd_box, 'r*', 'MarkerFaceColor','red', 'LineWidth', 1);
    hold on
    h2=errorbar(fres, dK_gss, dK_sd_gss, 'b^', 'MarkerFaceColor','blue', 'LineWidth', 1);
    
    legend([h1,h2],{'Box', 'Gaussian'},'interpreter','none','location','southeast')
%     legend({'unfiltered', ...
%         'box-filtered', ...
%         'Gaussian-filtered'})
        xlim([2 28])
        ylim([-0.35 0.05]);
    xlabel('\kappa', 'FontName', font, 'FontSize', 1.25*fontSize,'interpreter','tex')
    ylabel('$\frac{\delta (KE)}{KE}$','FontName',font,'FontSize',1.5*fontSize)
    title('(b) $\frac{\delta u}{u_0}$ = 1.5','FontName', font, 'FontSize', fontSize,'interpreter','latex','fontweight','normal')
end

if ismember('mag', display_plots)
    % Mean error plot.
%     figure;
%     errorbar(fres, dKm, dKm_sd, 'ko', 'MarkerFaceColor','black', 'LineWidth', 1)
%     hold on
    errorbar(fres, dKm_box, dKm_sd_box, 'r*', 'MarkerFaceColor','red', 'LineWidth', 1)
    hold on
    errorbar(fres, dKm_gss, dKm_sd_gss, 'b^', 'MarkerFaceColor','blue', 'LineWidth', 1)
    xlim([2 28])
    ylim([-0.025 0.35]);
    legend({'Box', 'Gaussian'},'interpreter','none')
%     legend({'unfiltered', ...
%         'box-filtered', ...
%         'Gaussian-filtered'})
    xlabel('\kappa', 'FontName', font, 'FontSize', 1.25*fontSize,'interpreter','tex')
    ylabel('$\frac{|\delta (KE)|}{KE}$','FontName',font,'FontSize',1.5*fontSize)
    title('(b) $\frac{\delta u}{u_0}$ = 1.5','FontName', font, 'FontSize', fontSize,'interpreter','latex','fontweight','normal')
end
