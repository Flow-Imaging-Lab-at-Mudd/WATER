function [err_mean_box, err_mean_gss, smoother_bias_box, smoother_bias_gss, ...
    amp_box, amp_gss] = KE_mean_err_run(vf, props)

range = vf.getRange();

% Each velocity component associated with a unit cell.
vol = prod(range(:,2) - range(:,1) + 1)*vf.solver.dv;
    
% Max speed in region of interest.
speed_r = vf.subsetVector(vf.data.speed);
max_speed_r = max(speed_r, [], 'all');
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

% Plot energy estimation error for small and large values of noise.
for i = 1: size(props, 2)
    vf.clearNoise();
    N = vf.noise_uniform(props(i)*u_mean);
    dK(i) = vf.kineticEnergy(1) - k;
    % Result with box smoothing.
    vf.smoothNoise('box');
    dK_box(i) = vf.kineticEnergy(1) - k;
    % Reset and smooth with gaussian filter.
    vf.setNoise(N)
    vf.smoothNoise('gaussian');
    dK_gss(i) = vf.kineticEnergy(1) - k;
end

% Normalize.
dK = dK / k;
dK_box = dK_box / k;
dK_gss = dK_gss / k;

% Baseline smoother biases. Assume 'prop' includes 0.
smoother_bias_box = abs(dK_box(1));
% Smoothing error on original velocity field by Gaussian.
smoother_bias_gss = abs(dK_gss(1));


% Mean smoothing errors.
err_mean_box = mean(abs(dK_box));

err_mean_gss = mean(abs(dK_gss));



% Smoothing errors as proportion of smoother bias.
err_prop_box = abs(dK_box / smoother_bias_box);
err_prop_gss = abs(dK_gss / smoother_bias_gss);
% Fit and record quadratic curves. The quadratic coefficient can serve as a
% measure of the KE error amplification rate.
err_quad_box = polyfit(props, err_prop_box, 2);
err_quad_gss = polyfit(props, err_prop_gss, 2);

% Amplification ratios.
amp_box = err_quad_box(1);
amp_gss = err_quad_gss(1);