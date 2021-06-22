function [dK, dKd, bias, err_mean] = ...
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
% Error from downsampling.
dKd = zeros(size(props));

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
end

% Normalize by correct KE.
dK = dK / k0;
dKd = dKd / k0;

abs_dK = abs(dK);
abs_dKd = abs(dKd);

% Baseline smoother biases.
bias = dKd(1);

abs_bias = abs(bias);

% Mean error of downsampling.
err_mean = mean(dKd);

abs_err_mean = abs(err_mean);

%%%%%%%%%%%%%%%%% Plot KE error %%%%%%%%%%%%%%%%%%
% figure;
% scatter(props, dK, 'filled')
% hold on
% scatter(props, dKd, 'r', 'filled')
% hold on
% yline(bias, '-')
% hold on
% yline(err_mean, '-')
% 
% 
% legend('unfiltered', 'downsampled', ...
%     strcat('downsampling bias $\kappa = $', {' '}, string(bias)),  ...
%     strcat('mean downsampling $\frac{\delta K}{K} = $', string(err_mean)), ...
%     'Interpreter', 'latex')
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\delta K}{K}$')

%%%%%%%%%%%%%%%%%% Plot absolute KE error %%%%%%%%%%%%%%%%%%%%%
figure;
scatter(props, abs_dK, 'filled')
hold on
scatter(props, abs_dKd, 'r', 'filled')
hold on
yline(abs_bias, '-')
hold on
yline(abs_err_mean, '-')

legend('unfiltered', 'downsampled', ...
    strcat('downsampling bias $\kappa = $', {' '}, string(abs(bias))),  ...
    strcat('mean downsampling $\left|\frac{\delta K}{K}\right| = $', string(abs(err_mean))), ...
    'Interpreter', 'latex')
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('$|\frac{\delta K}{K}|$')

%%%%%%%%%%%% KE error plot vs KE noise. %%%%%%%%%%%%
figure;
scatter(abs_dK, abs_dKd, 'r', 'filled')
hold on
% 1-1 line.
plot(abs_dK, abs_dK, 'black')
hold on
yline(abs_bias, '-', 'Color', 'r')

legend('downsampled', ...
    'identity line', ...
    strcat('downsampling bias $\kappa = $', string(abs_bias)), ...
    'Interpreter', 'latex')

xlabel('Unfiltered $|\frac{\delta K}{K}|$')
ylabel('Downsampled $\frac{|\delta K|}{\bar{K}}$')

%%%%%%%%%%%%% Smoothing errors as proportion of smoother bias %%%%%%%%%%%%
err_prop = abs(dKd / bias);
% Fit and record quadratic curves. The quadratic coefficient can serve as a
% measure of the KE error amplification rate.
err_quad = polyfit(props, err_prop, 2);

figure;
scatter(props, err_prop, 'r', 'filled')
hold on
err_fit = polyplot(err_quad, props);


legend(strcat('downsampling bias $\kappa = $', string(abs_bias)), ...
    strcat('fit $r^2 = $', string(cor(err_fit, err_prop))), ...
    'Interpreter', 'latex')
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
