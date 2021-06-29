function [dI, dI_box, dI_gss, bias_box, bias_gss] = ...
    impulse_err_run(vf, props, origin, fr, u0)

% Introduce levels of noise proportional to the mean speed in the effective
% region, according to 'props', e.g. 0: 0.1: 3. 'vf' is presume to have
% range properly set. 'origin' specifies the reference point to which
% impulse calculations are performed. 'u0' provides the freestream speed by
% which the theoretical impulse is obtained.
%
% Derek Li, June 2021

range = vf.getRange();
% Formatted string for title.
range_str = strcat('Range:', {' '}, mat2str(range));

props_count = length(props);

% Each velocity component associated with a unit cell.
vol = prod(range(:,2) - range(:,1) + 1)*vf.solver.dv;

% Set constant maximal magnitude of noise.
u_mean = vf.meanSpeed(0, 0);

% Theoretical momentum.
I0 = vf.fluid.density*[0 2*pi*fr^3*u0*vf.scale.len^4 0]';
i0 = norm(I0);
% I0 = vf.impulse(0, origin);
% i0 = norm(I0);

% Error in impulse computation given noise.
dI = zeros(3, props_count);
% Box smoothing.
dI_box = zeros(3, props_count);
% Gaussian smoothing.
dI_gss = zeros(3, props_count);

% Plot energy estimation error for small and large values of noise.
for i = 1: props_count
    vf.clearNoise();
    N = vf.noise_uniform(props(i)*u_mean);
    dI(:, i) = vf.impulse(1, origin) - I0;
    % Result with box smoothing.
    vf.smoothNoise('box');
    dI_box(:, i) = vf.impulse(1, origin) - I0;
    % Reset and smooth with gaussian filter.
    vf.setNoise(N)
    vf.smoothNoise('gaussian');
    dI_gss(:, i) = vf.impulse(1, origin) - I0;
end

% Normalize by magnitude of impulse in the region.
dI = dI / i0;
di = sqrt(sum(dI.^2, 1));
dI_box = dI_box / i0;
di_box = sqrt(sum(dI_box.^2, 1));
dI_gss = dI_gss / i0;
di_gss = sqrt(sum(dI_gss.^2, 1));

abs_dI = abs(dI);
abs_dI_box = abs(dI_box);
abs_dI_gss = abs(dI_gss);

% Baseline smoother biases.
bias_box = dI_box(:, 1);
bias_gss = dI_gss(:, 1);

mag_bias_box = norm(bias_box);
mag_bias_gss = norm(bias_gss);

abs_bias_box = abs(bias_box);
abs_bias_gss = abs(bias_gss);

%%%%%%%%%%% Plot signed impulse error %%%%%%%%%%%%%%%
% % y dimension.
% figure;
% scatter(props, dI(2,:))
% hold on
% scatter(props, dI_box(2,:), 'r', 'filled')
% hold on
% err_mean_box = mean(dI_box(2,:));
% yline(err_mean_box, '-')
% hold on
% scatter(props, dI_gss(2,:), 'b', 'filled')
% hold on
% err_mean_gss = mean(dI_gss(2,:));
% yline(err_mean_gss, '-')
% 
% legend({'unfiltered error', ...
%     'box-filtered $\vec{u}$', ...
%     strcat('box mean $\frac{\delta I_y}{I} = $', string(err_mean_box)), ...
%     'Gaussian-filtered $\vec{u}$', ...
%     strcat('Gaussian mean $\frac{\delta I_y}{I} = $', string(err_mean_gss))})
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\delta I_y}{I}$')
% title('Axial Impulse Error')
% 
% 
% % z dimension.
% figure;
% scatter(props, dI(3,:))
% hold on
% scatter(props, dI_box(3,:), 'r', 'filled')
% hold on
% err_mean_box = mean(dI_box(3,:));
% yline(err_mean_box, '-')
% hold on
% scatter(props, dI_gss(3,:), 'b', 'filled')
% hold on
% err_mean_gss = mean(dI_gss(3,:));
% yline(err_mean_gss, '-')
% 
% legend({'unfiltered error', ...
%     'box-filtered $\vec{u}$', ...
%     strcat('box mean $\frac{\delta I_z}{I} = $', string(err_mean_box)), ...
%     'Gaussian-filtered $\vec{u}$', ...
%     strcat('Gaussian mean $\frac{\delta I_z}{I} = $', string(err_mean_gss))})
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\delta I_z}{I}$')
% title('Impulse Error in $z$ Direction')
% 
% % x dimension.
% figure;
% scatter(props, dI(1,:))
% hold on
% scatter(props, dI_box(1,:), 'r', 'filled')
% hold on
% err_mean_box = mean(dI_box(1,:));
% yline(err_mean_box, '-')
% hold on
% scatter(props, dI_gss(1,:), 'b', 'filled')
% hold on
% err_mean_gss = mean(dI_gss(1,:));
% yline(err_mean_gss, '-')
% 
% legend({'unfiltered error', ...
%     'box-filtered $\vec{u}$', ...
%     strcat('box mean $\frac{\delta I_x}{I} = $', string(err_mean_box)), ...
%     'Gaussian-filtered $\vec{u}$', ...
%     strcat('Gaussian mean $\frac{\delta I_x}{I} = $', string(err_mean_gss))})
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\delta I_x}{I}$')
% title('Impulse Error in $x$ Direction')


%%%%%%%%%%%%%%%%%% Plot absolute impulse error %%%%%%%%%%%%%%%%%%%%

% % y dimension.
% figure;
% scatter(props, abs_dI(2,:))
% hold on
% err_mean0 = mean(abs_dI(2, :));
% yline(err_mean0, '-')
% hold on
% scatter(props, abs_dI_box(2,:), 'r', 'filled')
% hold on
% err_mean_box = mean(abs_dI_box(2,:));
% yline(err_mean_box, '-')
% hold on
% scatter(props, abs_dI_gss(2,:), 'b', 'filled')
% hold on
% err_mean_gss = mean(abs_dI_gss(2,:));
% yline(err_mean_gss, '-')
% 
% legend({'unfiltered error', ...
%     strcat('unfiltered mean $|\frac{\delta I_y}{I}| = $', string(err_mean0)), ...
%     'box-filtered $\vec{u}$', ...
%     strcat('box mean $|\frac{\delta I_y}{I}| = $', string(err_mean_box)), ...
%     'Gaussian-filtered $\vec{u}$', ...
%     strcat('Gaussian mean $|\frac{\delta I_y}{I}| = $', string(err_mean_gss))})
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\left|\frac{\delta I_y}{I}\right|$')
% title('Absolute Axial Impulse Error')
% 
% % z dimension.
% figure;
% scatter(props, abs_dI(3,:))
% hold on
% err_mean0 = mean(abs_dI(3, :));
% yline(err_mean0, '-')
% hold on
% scatter(props, abs_dI_box(3,:), 'r', 'filled')
% hold on
% err_mean_box = mean(abs_dI_box(3,:));
% yline(err_mean_box, '-')
% hold on
% scatter(props, abs_dI_gss(3,:), 'b', 'filled')
% hold on
% err_mean_gss = mean(abs_dI_gss(3,:));
% yline(err_mean_gss, '-')
% 
% legend({'unfiltered error', ...
%     strcat('unfiltered mean $|\frac{\delta I_z}{I}| = $', string(err_mean0)), ...
%     'box-filtered $\vec{u}$', ...
%     strcat('box mean $|\frac{\delta I_z}{I}| = $', string(err_mean_box)), ...
%     'Gaussian-filtered $\vec{u}$', ...
%     strcat('Gaussian mean $|\frac{\delta I_z}{I}| = $', string(err_mean_gss))})
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\left|\frac{\delta I_z}{I}\right|$')
% title('Absolute Impulse Error in $z$ Direction')
% 
% % x dimension.
% figure;
% scatter(props, abs_dI(1,:))
% hold on
% err_mean0 = mean(abs_dI(1, :));
% yline(err_mean0, '-')
% hold on
% scatter(props, abs_dI_box(1,:), 'r', 'filled')
% hold on
% err_mean_box = mean(abs_dI_box(1,:));
% yline(err_mean_box, '-')
% hold on
% scatter(props, abs_dI_gss(1,:), 'b', 'filled')
% hold on
% err_mean_gss = mean(abs_dI_gss(1,:));
% yline(err_mean_gss, '-')
% 
% legend({'unfiltered error', ...
%     strcat('unfiltered mean $|\frac{\delta I_x}{I}| = $', string(err_mean0)), ...
%     'box-filtered $\vec{u}$', ...
%     strcat('box mean $|\frac{\delta I_x}{I}| = $', string(err_mean_box)), ...
%     'Gaussian-filtered $\vec{u}$', ...
%     strcat('Gaussian mean $|\frac{\delta I_x}{I}| = $', string(err_mean_gss))})
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\left|\frac{\delta I_x}{I}\right|$')
% title('Absolute Impulse Error in $x$ Direction')

%%%%%%%%%%% Error magnitude vs proportional noise %%%%%%%%%%%%
% figure;
% scatter(props, di)
% hold on
% scatter(props, di_box, 'r', 'filled')
% hold on
% scatter(props, di_gss, 'b', 'filled')
% hold on
% yline(mag_bias_box, '-', 'Color', 'r')
% hold on
% yline(mag_bias_gss, '-', 'Color', 'b')
% 
% legend('unfiltered error', 'box-filtered', 'Gaussian-filtered', ...
%     strcat('box bias $\kappa = $', string(mag_bias_box)), ...
%     strcat('Gaussian bias $\kappa = $', string(mag_bias_gss)), ...
%     'Interpreter', 'latex')
% 
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('Filtered $\frac{|\delta I|}{\bar{I}}$')
% title('Smoother Efficacy')


%%%%%%%%%%%% Impulse error vs impulse noise, in magnitude%%%%%%%%%%%%
% figure;
% scatter(di, di_box, 'r', 'filled')
% hold on
% scatter(di, di_gss, 'b', 'filled')
% hold on
% % 1-1 line.
% plot(di, di, 'black')
% hold on
% yline(mag_bias_box, '-', 'Color', 'r')
% hold on
% yline(mag_bias_gss, '-', 'Color', 'b')
% 
% legend('box-filtered', 'Gaussian-filtered', ...
%     'identity line', ...
%     strcat('box bias $\kappa = $', string(mag_bias_box)), ...
%     strcat('Gaussian bias $\kappa = $', string(mag_bias_gss)), ...
%     'Interpreter', 'latex')
% 
% xlabel('Unfiltered $|\frac{\delta I}{I}|$')
% ylabel('Filtered $\frac{|\delta I|}{\bar{I}}$')
% title('Smoother Efficacy')

%%%%%%%%%%%%% Smoothing errors as proportion of smoother bias %%%%%%%%%%%%
% err_prop_box = abs(dI_box / mag_bias_box);
% err_prop_gss = abs(dI_gss / mag_bias_gss);
% 
% % y dimension.
% figure;
% scatter(props, err_prop_box(2,:), 'r', 'filled')
% hold on
% scatter(props, err_prop_gss(2,:), 'b', 'filled')
% 
% legend({strcat('box $\kappa = $', string(mag_bias_box)), ...
%     strcat('Gaussian $\kappa = $', string(mag_bias_gss))}, ...
%     'Interpreter', 'latex')
% 
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\left|\frac{\delta I}{I}\right|}{\kappa}$')
% title(strcat('Axial Error Proportional to Smoother Bias'))
% 
% % z dimension.
% figure;
% scatter(props, err_prop_box(3,:), 'r', 'filled')
% hold on
% scatter(props, err_prop_gss(3,:), 'b', 'filled')
% 
% legend({strcat('box $\kappa = $', string(mag_bias_box)), ...
%     strcat('Gaussian $\kappa = $', string(mag_bias_gss))}, ...
%     'Interpreter', 'latex')
% 
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\left|\frac{\delta I}{I}\right|}{\kappa}$')
% title(strcat('$z$ Error Proportional to Smoother Bias'))
% 
% 
% % x dimension.
% figure;
% scatter(props, err_prop_box(1,:), 'r', 'filled')
% hold on
% scatter(props, err_prop_gss(1,:), 'b', 'filled')
% 
% legend({strcat('box $\kappa = $', string(mag_bias_box)), ...
%     strcat('Gaussian $\kappa = $', string(mag_bias_gss))}, ...
%     'Interpreter', 'latex')
% 
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\left|\frac{\delta I}{I}\right|}{\kappa}$')
% title(strcat('$x$ Error Proportional to Smoother Bias'))
% 
