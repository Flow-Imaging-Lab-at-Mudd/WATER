vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
speed = sqrt(sum(vf.U.^2, 4));

% Minimal and maximal volume dimensions.
vol = [8 10; 8 10; 8 10];

% Introduce noise proportionally.
props = (1:30) * 0.1;

% Randomly sample effective regions.
num_ite = 20;
% Linear correlations.
cor = zeros([num_ite 2]);

for j = 1: num_ite
    % Random range.
    range = randRange(vf.dims, vol);
    range
    vf.setRange(range)
    
    % Max speed in region of interest.
    speed_r = vf.subsetVector(speed);
    max_speed_r = max(speed_r, [], 'all');
    % Set constant maximal magnitude of noise.
    err_mag = vf.meanSpeed(0) / 4;
    dK = zeros(size(props));
    
    % Kinetic energy without noise.
    k = vf.kineticEnergy(0);
    % Plot energy estimation error for small and large values of noise.
    for i = 1: size(props, 2)
        vf.clearNoise();
        vf.noise_uniform(props(i)*err_mag);
        dK(i) = vf.kineticEnergy(1) - k;
    end
    
    % Plot KE error.
    figure;
    subplot(2, 2, 1)
    scatter(props(1:10), dK(1:10), 'filled')
    xlabel('$|\delta u|$ (per mean speed)')
    ylabel('$\delta K$ (J)')
    subplot(2, 2, 2)
    scatter(props, dK, 'filled')
    xlabel('$|\delta u|$ (per mean speed)')
    ylabel('$\delta K$ (J)')
    
    % Plot absolute KE error.
    subplot(2, 2, 3)
    scatter(props(1:10), abs(dK(1:10)), 'filled')
    cor(i, 1) = corr(props(1:10), abs(dK(1:10)));
    title(strcat('$r = $', string(cor(i, 1))))
    xlabel('$|\delta u|$ (per mean speed)')
    ylabel('$|\delta K|$ (J)')
    subplot(2, 2, 4)
    scatter(props, abs(dK), 'filled')
    cor(i, 2) = corr(props, abs(dK));
    title(strcat('$r = $', string(cor(i, 2))))
    xlabel('$|\delta u|$ (per mean speed)')
    ylabel('$|\delta K|$ (J)')
    
    sgtitle(strcat('Volume:', {' '}, strjoin(string(range(:,2)-range(:,1))), ...
        {'; '}, '$\bar{u} = $', string(vf.meanSpeed(0))))
    
    uerr_histogram(vf.N_e);
    
%     vf.plotVector(vf.U_e, 0, '')
%     vf.plotScalar(sqrt(sum(vf.N_e.^2, 4)), 0, '')
    
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
