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
    
    % Run script for KE error samrpling.
    KE_uniform_err_run
    
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
