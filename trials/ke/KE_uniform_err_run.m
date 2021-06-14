% Presumed parameters: 'range', 'vf' with range properly set,
% 'vf.data.speed' for global speed without noise.

% Introduce noise proportionally.
props = 0: 0.1: 3;
% Index by which linearity no longer approximates.
lin_index = 5;

% Each velocity component associated with a unit cell.
vol = prod(range(:,2) - range(:,1) + 1)*vf.solver.dv;
    
% Max speed in region of interest.
speed_r = vf.subsetVector(vf.data.speed);
max_speed_r = max(speed_r, [], 'all');
% Set constant maximal magnitude of noise.
err_mag = vf.meanSpeed(0, 0);

% Kinetic energy without noise.
k = vf.kineticEnergy(0);

% Error in energy estimation given noise.
dK = zeros(size(props));
% Smoothed velocity results.
dK_us = zeros(size(props));

% Plot energy estimation error for small and large values of noise.
for i = 1: size(props, 2)
    vf.clearNoise();
    vf.noise_uniform(props(i)*err_mag);
    dK(i) = vf.kineticEnergy(1) - k;
    % Result with smoothing.
    vf.N_e(:,:,:,1) = smooth3(vf.U_e(:,:,:,1) + vf.N_e(:,:,:,1), 'box') - vf.U_e(:,:,:,1);
    vf.N_e(:,:,:,2) = smooth3(vf.U_e(:,:,:,2) + vf.N_e(:,:,:,2), 'box') - vf.U_e(:,:,:,2);
    vf.N_e(:,:,:,3) = smooth3(vf.U_e(:,:,:,3) + vf.N_e(:,:,:,3), 'box') - vf.U_e(:,:,:,3);
    dK_us(i) = vf.kineticEnergy(1) - k;
end

% Normalize.
dK = dK / k;
dK_us = dK_us / k;

% Formatted string for title.
range_str = strcat('Range:', {' '}, mat2str(range));

% % Plot KE error.
% figure;
% scatter(props, dK, 'filled')
% hold on
% scatter(props, dK_us, 'r', 'filled')
% hold on
% err_mean = mean(dK_us);
% yline(err_mean, '-')
% legend('unfiltered error', 'filtered $\vec{u}$', ...
%     strcat('$\frac{\delta K}{dK} = $', string(err_mean)), 'Interpreter', 'latex')
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\delta K}{K}$')
% title(range_str)

% Plot absolute KE error.
figure;
scatter(props, abs(dK), 'filled')
hold on
scatter(props, abs(dK_us), 'r', 'filled')
hold on
abs_err_mean = mean(abs(dK_us));
yline(abs_err_mean, '-')
legend('unfiltered error', 'filtered $\vec{u}$', ...
    strcat('$\frac{\delta K}{dK} = $', string(abs_err_mean)), 'Interpreter', 'latex')
title(range_str)

% Theoretical quadratic correlation.
pred = vf.fluid.density*vol*err_mag^2*vf.scale.len^2*(props + 1/2*props.^2) / k;
% plot(props, pred)
% title(strcat('$r = $', string(cor(i, 2))))

xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('$|\frac{\delta K}{K}|$')

% uerr_histogram(vf.N_e);

vf.plotVector(vf.U_e, 0, range_str);

%     vf.plotScalar(sqrt(sum(vf.N_e.^2, 4)), 0, '');
% plane_range = range;
% % Plot a parallel xy plane.f
% plane_range(3, 2) = plane_range(3, 1);
% vf.plotPlaneScalar(sqrt(sum(vf.N.^2, 4)), plane_range, 0, 'noise $\Delta u$')