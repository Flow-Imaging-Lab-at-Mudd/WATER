function [dK, dKd, dK_box, dK_gss, bias_box, bias_gss, abs_err_mean_box, abs_err_mean_gss] = ...
    KE_err_window_run(vf0, winsize, overlap, props, fr)
% Identical to KE_err_run.m except that error now by comparison of the
% downsampled data to the correct original data 'vf0'. The accurate
% velocity field is first downsampled, and then noise is introduced to the
% downsampled field.


% Presumed parameters: 'range', 'vf' with range properly set; and
% occasionally 'fr' when a mean central speed is to be estimated for a
% Hill's vortex.

% Correct kinetic energy computed from original data.
k0 = vf0.kineticEnergy(0);

% Create downsampled field and compute error.
vfd = vf0.downsample(winsize, overlap, 0);
% Shift between noise-free KE vals.
kd = vfd.kineticEnergy(0);
% KE error from downsampling.
dKd = kd - k0;

% Introduce error to the downsampled field. The positions are not altered
% in downsampling.
[dK, dK_box, dK_gss, ~, ~] = KE_err_run(vfd, props, fr);

% Refashion errors with respect to KE from original field.
dK = (kd*dK + dKd) / k0;
abs_dK = abs(dK);

dK_box = (kd*dK_box + dKd) / k0;
abs_dK_box = abs(dK_box);

dK_gss = (kd*dK_gss + dKd) / k0;
abs_dK_gss = abs(dK_gss);

bias_box = dK_box(1);
bias_gss = dK_gss(1);

err_mean_box = mean(dK_box);
abs_err_mean_box = mean(abs(dK_box));

err_mean_gss = mean(dK_gss);
abs_err_mean_gss = mean(abs(dK_gss));

% Now normalize sampling bias.
dKd = dKd / k0;


%%%%%%%%%%%%%%%%% Plot KE error %%%%%%%%%%%%%%%%%%
% figure;
% scatter(props, dK, 'filled')
% hold on
% yline(dKd, 'g', '-')
% hold on
% scatter(props, dK_box, 'r', 'filled')
% hold on
% scatter(props, dK_gss, 'b', 'filled')
% 
% legend('unfiltered downsampled error', ...
%     strcat('downsampling bias $\delta = $', {' '}, string(dKd)), ...
%     strcat('box-filtered, mean $\frac{\delta K}{K} = $', {' '}, string(err_mean_box)),  ...
%     strcat('box-filtered, mean $\frac{\delta K}{K} = $', {' '}, string(err_mean_gss)), ...
%     'Interpreter', 'latex')
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\delta K}{K}$')
% title(strcat('KE Error with Downsampling $\delta = $', ...
%     {' '}, string(abs(dKd))))

%%%%%%%%%%%%%%%%%% Plot absolute KE error %%%%%%%%%%%%%%%%%%%%%
% figure;
% scatter(props, abs(dK), 'filled')
% hold on
% yline(abs(dKd), 'g', '-')
% hold on
% scatter(props, abs(dK_box), 'r', 'filled')
% hold on
% scatter(props, abs(dK_gss), 'b', 'filled')
% 
% legend('unfiltered downsampled error', ...
%     strcat('downsampling bias $\delta = $', {' '}, string(abs(dKd))), ...
%     strcat('box-filtered, mean $\left|\frac{\delta K}{K}\right| = $', {' '}, string(abs_err_mean_box)),  ...
%     strcat('box-filtered, mean $\left|\frac{\delta K}{K}\right| = $', {' '}, string(abs_err_mean_gss)), ...
%     'Interpreter', 'latex')
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\left|\frac{\delta K}{K}\right|$')
% title(strcat('Absolute KE Error with Downsampling $\delta = $', ...
%     {' '}, string(abs(dKd))))
% 
%%%%%%%%%%%%%% KE error plot vs KE noise. %%%%%%%%%%%%
% figure;
% yline(abs(dKd), '-', 'Color', 'g')
% hold on
% % 1-1 line.
% plot(abs_dK, abs_dK, 'black')
% hold on
% scatter(abs_dK, abs_dK_box, 'r', 'filled')
% hold on
% scatter(abs_dK, abs_dK_gss, 'b', 'filled')
% 
% legend('downsampling bias $\delta$', ...
%     'identity line', ...
%     strcat('box-filtered, bias $\kappa = $', string(abs(bias_box))), ...
%     strcat('Gaussian-filtered, bias $\kappa = $', string(abs(bias_gss))), ...
%     'Interpreter', 'latex')
% 
% xlabel('Unfiltered $|\frac{\delta K}{K}|$')
% ylabel('$|\frac{\delta K}{\bar{K}}|$')
% title(strcat('Filerting of KE Noise with Downsampling $\delta = $', ...
%     {' '}, string(abs(dKd))))

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
% hold on
% scatter(props, err_prop_gss, 'b', 'filled')
% hold on
% err_fit_gss = polyplot(err_quad_gss, props, 'b');
% 
% 
% legend(strcat('box-filtered, $\kappa = $', string(abs(bias_box))), ...
%     strcat('box fit $r^2 = $', string(cor(err_fit_box, err_prop_box))), ...
%     strcat('Gaussian-filtered, $\kappa = $', string(abs(bias_gss))), ...
%     strcat('Gaussian fit $r^2 = $', string(cor(err_fit_gss, err_prop_gss))), ...
%     'Interpreter', 'latex')
% 
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\left|\frac{\delta K}{K}\right|}{\kappa}$')
% title(strcat('Error Proportional to Smoother Bias with Downsampling $\delta = $', ...
%     {' '}, string(abs(dKd))))