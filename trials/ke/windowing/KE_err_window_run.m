function [dK, dKd, dK_box, dK_gss, bias_box, bias_gss] = ...
    KE_err_window_run(vf0, winsize, overlap, props, fr)
% Identical to KE_err_run.m except that error now by comparison of the
% downsampled data to the correct original data 'vf0'. Further, as the
% downsampling implemented involved averaging of velocities within window,
% we apply no additional filter.


% Presumed parameters: 'range', 'vf' with range properly set; and
% occasionally 'fr' when a mean central speed is to be estimated for a
% Hill's vortex.

range0 = vf0.getRange();

% % Set constant maximal magnitude of noise.
% u_mean = vf0.meanSpeed(0, 0);

% Use central speed for feature focusing. Enough resolution is presumed for
% proper indexing.
center = [floor(vf0.getDims()/2) - floor(fr / vf0.xresol); ...
    floor(vf0.getDims()/2) + floor(fr / vf0.xresol)]';
center(:, 1) = max([center(:,1) ones(3, 1)], [], 2);
center(:, 2) = min([center(:,2) vf0.getDims()'], [], 2);
vf0.setRange(center)
u_mean = vf0.meanSpeed(0, 0);
vf0.setRange(range0)

% Correct kinetic energy computed from original data.
k0 = vf0.kineticEnergy(0);

% Error in energy estimation given noise.
dK = zeros(size(props));
% Error from downsampling without smoothing.
dKd = zeros(size(props));
% Box smoothing.
dK_box = zeros(size(props));
% Gaussian smoothing.
dK_gss = zeros(size(props));

% % KE error profiles per point.
% dK_pro_box = zeros(vf.span);
% dK_pro_gss = zeros(vf.span);

% Plot energy estimation error for small and large values of noise.
for i = 1: size(props, 2)
    vf0.clearNoise();
    N = vf0.noise_uniform(props(i)*u_mean);
    
    % Create downsampled field and compute error.
    vfd = vf0.downsample(winsize, overlap, 1);
    % Original KE noise.
    dK(i) = vf0.kineticEnergy(1) - k0;
    % Downsampled KE noise.
    dKd(i) = vfd.kineticEnergy(0) - k0;
    
    % Results with smoothing.
    vfd.smoothVelocity('box');
    dK_box(i) = vfd.kineticEnergy(0) - k0;
%    dK_pro_box = abs(KE_profile(vf, 1) - KE_profile(vf, 0)) / k;
    vfd.smoothVelocity('gaussian');
    dK_gss(i) = vfd.kineticEnergy(0) - k0;
%    dK_pro_gss = abs(KE_profile(vf, 1) - KE_profile(vf, 0)) / k;
    
end

% Normalize by correct KE.
dK = dK / k0;
dKd = dKd / k0;
dK_box = dK_box / k0;
dK_gss = dK_gss / k0;


% Baseline smoother biases.
bias_box = dK_box(1);
bias_gss = dK_gss(1);

%%%%%%%%%%%%%%%%% Plot KE error %%%%%%%%%%%%%%%%%%
figure;
scatter(props, dK, 'filled')
hold on
scatter(props, dKd, 'g', 'filled')
hold on
scatter(props, dK_box, 'r', 'filled')
hold on
err_mean_box = mean(dK_box);
yline(err_mean_box, '-')
hold on
scatter(props, dK_gss, 'b', 'filled')
hold on
err_mean_gss = mean(dK_gss);
yline(err_mean_gss, '-')

legend('unfiltered $\delta K$', 'unfiltered \& downsampled', ...
    'box-filtered $\vec{u}$',  ...
    strcat('box $\frac{\delta K}{K} = $', string(err_mean_box)), ...
    'Gaussian-filtered $\vec{u}$', ...
    strcat('Gaussian $\frac{\delta K}{K} = $', string(err_mean_gss)), ...
    'Interpreter', 'latex')
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('$\frac{\delta K}{K}$')
% title(range_str)

%%%%%%%%%%%%%%%%%% Plot absolute KE error %%%%%%%%%%%%%%%%%%%%%
figure;
scatter(props, abs(dK), 'filled')
hold on
scatter(props, abs(dKd), 'g', 'filled')
hold on
scatter(props, abs(dK_box), 'r', 'filled')
hold on
abs_err_mean_box = mean(abs(dK_box));
yline(abs_err_mean_box, '-')
hold on
scatter(props, abs(dK_gss), 'b', 'filled')
hold on
abs_err_mean_gss = mean(abs(dK_gss));
yline(abs_err_mean_gss, '-')

legend('unfiltered $|\delta K|$', ...
    'downsampled', ...
    'box-filtered $\vec{u}$', ...
    strcat('box $\left|\frac{\delta K}{K}\right| = $', ...
    string(abs_err_mean_box)), ...
    'Gaussian-filtered $\vec{u}$', ...
    strcat('Gaussian $\left|\frac{\delta K}{K}\right| = $', string(abs_err_mean_gss)), ...
    'Interpreter', 'latex')
%     strcat('box bias $\kappa = $', string(abs(bias_box))), ...
%     strcat('Gaussian bias $\kappa = $', string(abs(bias_gss))), ...
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('$|\frac{\delta K}{K}|$')
% title(range_str)

%%%%%%%%%%%% KE error plot vs KE noise. %%%%%%%%%%%%
% figure;
% scatter(abs(dK), abs(dK_box), 'r', 'filled')
% hold on
% scatter(abs(dK), abs(dK_gss), 'b', 'filled')
% hold on
% % 1-1 line.
% plot(dK, dK, 'black')
% hold on
% % Smoother biases.
% abs_bias_box = abs(bias_box);
% yline(abs_bias_box, '-', 'Color', 'r')
% 
% hold on
% abs_bias_gss = abs(bias_gss);
% yline(abs_bias_gss, '-', 'Color', 'b')
% 
% legend('box-filtered', 'Gaussian-filtered', ...
%     '$y=x$ identity line', ...
%     strcat('box bias $\kappa = $', string(abs_bias_box)), ...
%     strcat('Gaussian bias $\kappa = $', string(abs_bias_gss)), ...
%     'Interpreter', 'latex')
% 
% xlabel('Unfiltered $|\frac{\delta K}{K}|$')
% ylabel('Filtered $\frac{|\delta K|}{\bar{K}}$')
% title(range_str)

%%%%%%%%%%%%% Smoothing errors as proportion of smoother bias %%%%%%%%%%%%
% err_prop_box = abs(dK_box / bias_box);
% err_prop_gss = abs(dK_gss / bias_gss);
% % Fit and record quadratic curves. The quadratic coefficient can serve as a
% % measure of the KE error amplification rate.
% err_quad_box = polyfit(props, err_prop_box, 2);
% err_quad_gss = polyfit(props, err_prop_gss, 2);
% 
% figure;
% scatter(props, err_prop_box, 'r', 'filled')
% hold on
% err_fit_box = polyplot(err_quad_box, props);
% 
% hold on
% scatter(props, err_prop_gss, 'b', 'filled')
% hold on
% err_fit_gss = polyplot(err_quad_gss, props);
% 
% 
% legend(strcat('box $\kappa = $', string(abs(bias_box))), ...
%     strcat('box fit $r^2 = $', string(cor(err_fit_box, err_prop_box))), ...
%     strcat('Gaussian $\kappa = $', string(abs(bias_gss))), ...
%     strcat('Gaussian fit $r^2 = $', string(cor(err_fit_gss, err_prop_gss))), ...
%     'Interpreter', 'latex')
% 
% % Different test dataset run.
% % dataset_str = 'Turbulent Vortex Ring';
% dataset_str = 'Synthetic Hill Vortex';
% 
% 
% title(strcat('Error Proportional to Smoother Bias:', {' '}, dataset_str))
% 
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\left|\frac{\delta K}{K}\right|}{\kappa}$')


% uerr_histogram(vf.N_e);

% vf.plotVector(vf.U_e, 0, strcat(range_str, {' '}, '$\bar{u} = $', string(u_mean)));

% vf.plotScalar(sqrt(sum(vf.N_e.^2, 4)), 0, '');

% plane_range = range;
% % Plot a parallel xy plane.f
% plane_range(3, 2) = plane_range(3, 1);
% vf.plotPlaneScalar(sqrt(sum(vf.N.^2, 4)), plane_range, 0, 'noise $\Delta u$')
