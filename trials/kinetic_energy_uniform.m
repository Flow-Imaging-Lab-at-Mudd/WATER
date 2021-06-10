vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
vf.setRange([27 30; 20 22; 10 12]);
vf.data.speed = sqrt(sum(vf.U.^2, 4));

% Maximal volume dimensions.
vol_max = [5 5 5];

% Introduce noise proportionally.
props = (1:30) * 0.1;

% Randomly sample effective regions.
num_ite = 20;
for j = 1: num_ite
    % Random range.
    range = [randi(vf.dims(2) - vol_max(1)); ...
        randi(vf.dims(1) - vol_max(2)); ...
        randi(vf.dims(3) - vol_max(3))];
    range(:, 2) = range(:, 1) + [randi(vol_max(1)); randi(vol_max(2)); randi(vol_max(3))];
    disp(range)
    vf.setRange(range)
    
    % Max speed in region of interest.
    speed_r = vf.subsetVector(vf.data.speed);
    max_speed_r = max(speed_r, [], 'all');
    % Set constant maximal magnitude of noise.
    err_mag = max_speed_r / 5;
    dK = zeros(size(props));
    
    % Kinetic energy without noise.
    k = vf.kineticEnergy(0);
    % Plot energy estimation error for small and large values of noise.
    for i = 1: size(props, 2)
        vf.clearNoise();
        vf.noise_uniform(props(i)*err_mag);
        dK(i) = vf.kineticEnergy(1) - k;
    end
    figure;
    subplot(1, 2, 1)
    scatter(props(1:10), dK(1:10))
    xlabel('$\delta u$ (per maximum speed)')
    ylabel('$\delta K$ (J)')
    subplot(1, 2, 2)
    scatter(props, dK)
    xlabel('$\delta u$ (per maximum speed)')
    ylabel('$\delta K$ (J)')
    title(strcat('Volume ', string(strjoin(string(range(:,2)-range(:,1))))))
    
    pause
    close
end
