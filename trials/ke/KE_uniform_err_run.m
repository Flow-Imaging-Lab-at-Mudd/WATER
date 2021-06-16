function [dK, dK_box, dK_gss, bias_box, bias_gss] = ...
    KE_uniform_err_run(vf, range)
% Presumed parameters: 'range', 'vf' with range properly set,
% 'vf.data.speed' for global speed without noise.

% Introduce noise proportionally.
props = 0: 0.1: 3;

% Each velocity component associated with a unit cell.
vol = prod(range(:,2) - range(:,1) + 1)*vf.solver.dv;
    
% Set constant maximal magnitude of noise.
u_mean = vf.meanSpeed(0, 0);

% Kinetic energy without noise.
k = vf.kineticEnergy(0);

% Error in energy estimation given noise.
dK = zeros(size(props));
% Box smoothing.
dK_box = zeros(size(props));
% Gaussian smoothing.
dK_gss = zeros(size(props));

% % KE error profiles per point.
% dK_pro_box = zeros(vf.span);
% dK_pro_gss = zeros(vf.span);

% Plot energy estimation error for small and large values of noise.
for i = 1: size(props, 2)
    vf.clearNoise();
    N = vf.noise_uniform(props(i)*u_mean);
    dK(i) = vf.kineticEnergy(1) - k;
    % Result with box smoothing.
    vf.smoothNoise('box');
    dK_box(i) = vf.kineticEnergy(1) - k;
%     dK_pro_box = abs(KE_profile(vf, 1) - KE_profile(vf, 0)) / k;
    % Reset and smooth with gaussian filter.
    vf.setNoise(N)
    vf.smoothNoise('gaussian');
    dK_gss(i) = vf.kineticEnergy(1) - k;
%     dK_pro_gss = abs(KE_profile(vf, 1) - KE_profile(vf, 0)) / k;
%     vf.plotScalar(dK_pro_box, 0, 'box');
%     vf.plotScalar(dK_pro_gss, 0, 'Gaussian');
end

% Normalize.
dK = dK / k;
dK_box = dK_box / k;
dK_gss = dK_gss / k;

% Baseline smoother biases.
bias_box = dK_box(1);
bias_gss = dK_gss(1);

% Formatted string for title.
range_str = strcat('Range:', {' '}, mat2str(range));

% % Plot KE error.
% % figure;
% % scatter(props, dK, 'filled')
% % hold on
% % scatter(props, dK_box, 'r', 'filled')
% % hold on
% % err_mean = mean(dK_box);
% % yline(err_mean, '-')
% % legend('unfiltered error', 'filtered $\vec{u}$', ...
% %     strcat('$\frac{\delta K}{K} = $', string(err_mean)), 'Interpreter', 'latex')
% % xlabel('$\frac{|\delta u|}{\bar{u}}$')
% % ylabel('$\frac{\delta K}{K}$')
% % title(range_str)
% 
% % Plot absolute KE error.
% figure;
% scatter(props, abs(dK), 'filled')
% hold on
% scatter(props, abs(dK_box), 'r', 'filled')
% hold on
% abs_err_mean_box = mean(abs(dK_box));
% yline(abs_err_mean_box, '-')
% hold on
% scatter(props, abs(dK_gss), 'y', 'filled')
% hold on
% abs_err_mean_gss = mean(abs(dK_gss));
% yline(abs_err_mean_gss, '-')
% 
% legend('unfiltered error', ...
%     'box-filtered $\vec{u}$', ...
%     strcat('box $\left|\frac{\delta K}{K}\right| = $', ...
%     string(abs_err_mean_box)), ...
%     'Gaussian-filtered $\vec{u}$', ...
%     strcat('Gaussian $\left|\frac{\delta K}{K}\right| = $', string(abs_err_mean_gss)), ...
%     'Interpreter', 'latex')
% %     strcat('box bias $\kappa = $', string(abs(bias_box))), ...
% %     strcat('Gaussian bias $\kappa = $', string(abs(bias_gss))), ...
% % Theoretical quadratic correlation.
% pred = vf.fluid.density*vol*u_mean^2*vf.scale.len^2*(props + 1/2*props.^2) / k;
% % plot(props, pred)
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$|\frac{\delta K}{K}|$')
% title(range_str)

% % Error plot vs KE noise.
% figure;
% scatter(abs(dK), abs(dK_box), 'r', 'filled')
% hold on
% scatter(abs(dK), abs(dK_gss), 'y', 'filled')
% 
% legend('box-filtered', 'Gaussian-filtered')
% xlabel('Unfiltered $|\frac{\delta K}{K}|$')
% ylabel('Filtered $\frac{|\delta K|}{\bar{K}}$')
% title(range_str)
% 
% % Smoothing errors as proportion of smoother bias.
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
% err_box_str = polyplot(err_quad_box, props);
% 
% hold on
% scatter(props, err_prop_gss, 'y', 'filled')
% hold on
% err_gss_str = polyplot(err_quad_gss, props);
% 
% 
% legend(strcat('box $\kappa = $', string(abs(bias_box))), ...
%     strcat('box fit'), ...
%     strcat('Gaussian $\kappa = $', string(abs(bias_gss))), ...
%     strcat('Gaussian fit'), ...
%     'Interpreter', 'latex')
% title('Proportional Error to Smoother Bias')
% 
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\left|\frac{\delta K}{K}\right|}{\kappa}$')




% uerr_histogram(vf.N_e);

% vf.plotVector(vf.U_e, 0, strcat(range_str, {' '}, '$\bar{u} = $', string(u_mean)));

%     vf.plotScalar(sqrt(sum(vf.N_e.^2, 4)), 0, '');
% plane_range = range;
% % Plot a parallel xy plane.f
% plane_range(3, 2) = plane_range(3, 1);
% vf.plotPlaneScalar(sqrt(sum(vf.N.^2, 4)), plane_range, 0, 'noise $\Delta u$')
