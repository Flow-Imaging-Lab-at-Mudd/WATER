vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
vf.data.speed = sqrt(sum(vf.U.^2, 4));

% High speed region.
% range_hs = [26 39; 7 22; 17 26];
% High-speed center.
range_hs = [31 36; 13 17; 20 24];

% Incident on one side. Differ only in z range.
% % Range enveloping the high-speed region. 
% range_en = [25 46; 9 29; 9 33];
% % Range in which the high-speed region is on the boundary.
% range_bound = [25 46; 9 29; 16 40];
% % Note that both regions enclose the high-speed region, which is completely
% % smoothed by both.

% Incident on two sides.
% Range enveloping the high-speed region. 
range_en = [25 43; 5 23; 9 33];
% Range in which the high-speed region is on the boundary.
range_bound = [25 43; 7 26; 16 40];
% Note that both regions enclose the high-speed region, which is completely
% smoothed by both.

% Overall region containing the three former regions on which noise is
% introduced.
range = [20 52; 5 32; 7 40];
vf.setRange(range)

%%%%%%%% Now introduce noise and consider smoothing of the high speed
%%%%%%%% region.

% Introduce noise proportionally.
props = 0: 0.1: 4;

% Each velocity component associated with a unit cell.
vol = prod(range(:,2) - range(:,1) + 1)*vf.solver.dv;
    
% Max speed in the overall region.
speed_r = vf.subsetVector(vf.data.speed);
max_speed_r = max(speed_r, [], 'all');
% Set constant maximal magnitude of noise.
u_mean_r = vf.meanSpeed(0, 0);

% Kinetic energy without noise in the high-speed region.
vf.setRange(range_hs)
k = vf.kineticEnergy(0);

% Error in energy estimation given noise.
dK = zeros(size(props));
% Error from smoothed velocity results from the enveloping region.
dK_en = zeros(size(props));
% Error from smoothed velocity results from the bounding region.
dK_bound = zeros(size(props));

% Plot energy estimation error for small and large values of noise.
for i = 1: size(props, 2)
    vf.clearNoise();
    vf.setRange(range)
    N_e = vf.noise_uniform(props(i)*err_mag);
    % All errors are computed on the high-speed region.
    vf.setRange(range_hs)
    dK(i) = vf.kineticEnergy(1) - k;
    % Smoothing in enveloping region.
    vf.setRange(range_en)
    vf.N_e(:,:,:,1) = smooth3(vf.U_e(:,:,:,1) + vf.N_e(:,:,:,1), 'box') - vf.U_e(:,:,:,1);
    vf.N_e(:,:,:,2) = smooth3(vf.U_e(:,:,:,2) + vf.N_e(:,:,:,2), 'box') - vf.U_e(:,:,:,2);
    vf.N_e(:,:,:,3) = smooth3(vf.U_e(:,:,:,3) + vf.N_e(:,:,:,3), 'box') - vf.U_e(:,:,:,3);
    % Update noise globally.
    vf.setNoise(vf.N_e)
    % Compute dKE.
    vf.setRange(range_hs)
    dK_en(i) = vf.kineticEnergy(1) - k;
    % Reset noise.
    vf.clearNoise()
    vf.setRange(range)
    vf.setNoise(N_e)
    % Smoothing in bounding region.
    vf.setRange(range_bound)
    vf.N_e(:,:,:,1) = smooth3(vf.U_e(:,:,:,1) + vf.N_e(:,:,:,1), 'box') - vf.U_e(:,:,:,1);
    vf.N_e(:,:,:,2) = smooth3(vf.U_e(:,:,:,2) + vf.N_e(:,:,:,2), 'box') - vf.U_e(:,:,:,2);
    vf.N_e(:,:,:,3) = smooth3(vf.U_e(:,:,:,3) + vf.N_e(:,:,:,3), 'box') - vf.U_e(:,:,:,3);
    % Update noise globally.
    vf.setNoise(vf.N_e);
    % Compute dKE.
    vf.setRange(range_hs)
    dK_bound(i) = vf.kineticEnergy(1) - k;
end

% Normalize.
dK = dK / k;
dK_en = dK_en / k;
dK_bound = dK_bound / k;

% % Plot KE error.
% figure;
% scatter(props, dK)
% hold on
% scatter(props, dK_en, 'filled')
% hold on
% scatter(props, dK_bound, 80, 'black')
% legend('unfiltered error', 'enveloping filter', 'boundary filter')
% xlabel('$\frac{|\delta u|}{\bar{u}}$')
% ylabel('$\frac{\delta K}{K}$')
% title('Error for Enveloping and Boundary Filters')

% Plot absolute KE error.
figure;
scatter(props, abs(dK), 'filled')
hold on
scatter(props, abs(dK_en), 'filled')
hold on
scatter(props, abs(dK_bound), 80, 'black')
hold on
abs_err_mean = mean(abs(dK_en));
yline(abs_err_mean, '-')
legend('unfiltered error', 'enveloping filter', 'boundary filter', ...
    strcat('$\frac{\delta K}{K} = $', string(abs_err_mean)), 'Interpreter', 'latex')
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('$|\frac{\delta K}{K}|$')
title('Error Magnitude in Enveloping Region')
