function [dK, dK_box, dK_gss, bias_box, bias_gss] = ...
    KE_err_run(vf, props, fr)
% Presumed parameters: 'range', 'vf' with range properly set; and
% occasionally 'fr' when a mean central speed is to be estimated for a
% Hill's vortex.

range = vf.getRange();

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
%    dK_pro_box = abs(KE_profile(vf, 1) - KE_profile(vf, 0)) / k;
    % Reset and smooth with gaussian filter.
    vf.setNoise(N)
    vf.smoothNoise('gaussian');
    dK_gss(i) = vf.kineticEnergy(1) - k;
%    dK_pro_gss = abs(KE_profile(vf, 1) - KE_profile(vf, 0)) / k;
%     vf.plotScalar(dK_pro_box, 0, 'box');
%     vf.plotScalar(dK_pro_gss, 0, 'Gaussian');
%     vf.plotVector(vf.U_e, 0, '$\vec{u}$')
    
%     pause
%     close all
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

% Plot KE error.
% figure;
% scatter(props, dK, 'filled')
% hold on
% scatter(props, dK_box, 'r', 'filled')
% hold on
% err_mean = mean(dK_box);
% yline(err_mean, '-')
% legend('unfiltered error', 'filtered $\vec{u}$', ...
%     strcat('$\frac{\delta K}{K} = $', string(err_mean)), 'Interpreter', 'latex')
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\delta K}{K}$')
% title('KE Error')

%%%%%%%%%%%%%%%%%% Plot absolute KE error %%%%%%%%%%%%%%%%%%%%%
% figure;
% scatter(props, abs(dK), 'filled')
% hold on
% scatter(props, abs(dK_box), 'r', 'filled')
% hold on
% abs_err_mean_box = mean(abs(dK_box));
% yline(abs_err_mean_box, '-')
% hold on
% scatter(props, abs(dK_gss), 'b', 'filled')
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
% 
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$|\frac{\delta K}{K}|$')
% title('Absolute KE Error')

% %     strcat('box bias $\kappa = $', string(abs(bias_box))), ...
% %     strcat('Gaussian bias $\kappa = $', string(abs(bias_gss))), ...
% % Theoretical quadratic correlation.
% pred = vf.fluid.density*vol*u_mean^2*vf.scale.len^2*(props + 1/2*props.^2) / k;
% % plot(props, pred)


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
% title('Filtering KE Noise')

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
% err_fit_box = polyplot(err_quad_box, props, 'r');
% 
% hold on
% scatter(props, err_prop_gss, 'b', 'filled')
% hold on
% err_fit_gss = polyplot(err_quad_gss, props, 'b');
% 
% 
% legend(strcat('box $\kappa = $', string(abs(bias_box))), ...
%     strcat('box fit $r^2 = $', string(cor(err_fit_box, err_prop_box))), ...
%     strcat('Gaussian $\kappa = $', string(abs(bias_gss))), ...
%     strcat('Gaussian fit $r^2 = $', string(cor(err_fit_gss, err_prop_gss))), ...
%     'Interpreter', 'latex')
% 
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\left|\frac{\delta K}{K}\right|}{\kappa}$')
% title(strcat('Error Proportional to Smoother Bias'))
