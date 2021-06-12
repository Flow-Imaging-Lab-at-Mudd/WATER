vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
speed = sqrt(sum(vf.U.^2, 4));

% Minimal and maximal volume dimensions.
vol_range = [20 30; 20 30; 20 30];

% Introduce noise proportionally.
props = 0: 0.1: 3;
% Index by which linearity no longer approximates.
lin_index = 5;

% Randomly sample effective regions.
num_ite = 20;
% Linear correlations.
cor = zeros([num_ite 2]);

for j = 1: num_ite
    % Random range.
    range = randRange(vf.dims, vol_range);
    range
    vf.setRange(range)
    % A copy of the regional velocity field, used for resetting later.
    U_e = vf.U_e;
    % Each velocity component associated with a unit cell.
    vol = prod(range(:,2) - range(:,1) + 1)*vf.solver.dv;
    
    % Max speed in region of interest.
    speed_r = vf.subsetVector(speed);
    max_speed_r = max(speed_r, [], 'all');
    % Set constant maximal magnitude of noise.
    err_mag = vf.meanSpeed(0, 0);
    dK = zeros(size(props));
    % Smoothed velocity results.
    dK_us = zeros(size(props));
    
    
    % Kinetic energy without noise.
    k = vf.kineticEnergy(0);
    % Plot energy estimation error for small and large values of noise.
    for i = 1: size(props, 2)
        vf.clearNoise();
        vf.noise_uniform(props(i)*err_mag);
        dK(i) = vf.kineticEnergy(1) - k;
        % Result with smoothing.
        vf.U_e(:,:,:,1) = smooth3(vf.U_e(:,:,:,1) + vf.N_e(:,:,:,1), 'box');
        vf.U_e(:,:,:,2) = smooth3(vf.U_e(:,:,:,2) + vf.N_e(:,:,:,2), 'box');
        vf.U_e(:,:,:,3) = smooth3(vf.U_e(:,:,:,3) + vf.N_e(:,:,:,3), 'box');
        dK_us(i) = vf.kineticEnergy(0) - k;
        vf.U_e = U_e;
    end
    
    % Normalize.
    dK = dK / k;
    dK_us = dK_us / k;
    
    % Plot KE error.
    figure;
    subplot(2, 2, 1)
    scatter(props(1:lin_index), dK(1:lin_index), 'filled')
    hold on
    scatter(props(1:lin_index), dK_us(1:lin_index), 'r', 'filled')
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    ylabel('$\frac{\delta K}{K}$')
    
    subplot(2, 2, 2)
    scatter(props, dK, 'filled')
    hold on
    scatter(props, dK_us, 'r', 'filled')
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    ylabel('$\frac{\delta K}{K}$')
    
    % Plot absolute KE error.
    subplot(2, 2, 3)
    scatter(props(1:lin_index), abs(dK(1:lin_index)), 'filled')
    hold on
    scatter(props(1:lin_index), abs(dK_us(1:lin_index)), 'r', 'filled')
    cor(i, 1) = corr(props(1:10), abs(dK(1:10)));
    % title(strcat('$r = $', string(cor(i, 1))))
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    ylabel('$|\frac{\delta K}{K}|$')
    
    subplot(2, 2, 4)
    scatter(props, abs(dK), 'filled')
    hold on
    scatter(props, abs(dK_us), 'r', 'filled')
    title(strcat('mean error percentage = ', {' '}, string(mean(abs(dK_us)))))
    % Theoretical quadratic correlation.
    pred = vf.fluid.density*vol*err_mag^2*vf.scale.len^2*(props + 1/2*props.^2) / k;
    % plot(props, pred)
    cor(i, 2) = corr(props, abs(dK));
    % title(strcat('$r = $', string(cor(i, 2))))
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    ylabel('$|\frac{\delta K}{K}|$')
    
%     sgtitle(strcat('$\frac{\bar{u}V}{}=$', {' '}, ...
%         string(prod(range(:,2)-range(:,1))/err_mag), ...
%         {' $\textrm{mm}^{2}$s'}))

    sgtitle(strcat('Volume:', {' '}, strjoin(string(range(:,2)-range(:,1)))))
    
    % uerr_histogram(vf.N_e);
    
    vf.plotVector(vf.U_e, 0, '');
%     vf.plotScalar(sqrt(sum(vf.N_e.^2, 4)), 0, '');
    plane_range = range;
    % Plot a parallel xy plane.f
    plane_range(3, 2) = plane_range(3, 1);
    % vf.plotPlaneScalar(sqrt(sum(vf.N.^2, 4)), plane_range, 0, 'noise $\Delta u$')
    
    pause
    close all
end

% Short hand for computing correlations given two row vectors.
function cor = corr(r1, r2)
    cor = corrcoef([r1; r2]');
    cor = cor(2, 1);
end

function plt = uerr_histogram(N)
    plt = figure;
    N = sqrt(sum(N.^2, 4));
    histogram(N(:));
    xlabel('$\Delta u$')
    ylabel('frequency')
end
