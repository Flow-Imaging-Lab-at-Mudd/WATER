% Track objective origin over time frames and assess consistency of impulse
% evaluation.
%
% Derek Li, July 2021

% Load data from folder.
fdr = 'C:\Users\derek\flow\data\L18_Run1\';
data = dir(fdr);

% Subfolder names corresponding to frames.
frame_fdrs = strings(1, size(data, 1));
% Import subfolders in sorted order.
for i = 1: size(data, 1)
    frame_fdrs(i) = strcat(data(i).name, '\');
end
frame_fdrs = sort(frame_fdrs);
% First two entries not valid subfolders. Next two with null velocity
% fields.
frame_fdrs = frame_fdrs(5:end);

% Number of time frames to consider.
num_frame = 10;
% Array of velocity fields, corresonding to different time frames.
vfs = cell(1, num_frame);
% Objective origin obtained in volume restricted to the vortical region.
origin_obj_res = zeros(3, num_frame);
% Objective origin obtained in the global volume.
origin_obj_glo = zeros(3, num_frame);

% Track objective origin over time. The motion of the feature vortex is
% estimated empirically and hardcoded here. A plot produced here is used to
% help assess the fidelity.
for i = 1: num_frame
    % Advance by multiple frames to visualize shift.
    load(strcat(fdr, frame_fdrs(-2 + 3*i), '\3DPIV_postprocessed_results_calibrated.mat'), ...
        'x', 'y', 'z', 'u', 'v', 'w')
    vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
    % Tentative position of center of vortex. Determined manually. The 'x'
    % and 'z' positions of the vortex do not alter.
    center = [-8 10-1.5*(i-1) -21];
    % Enclose vortical region, whose dimension empirically does not change appreciably.
    % vf.setRangePosition([-20 0; center(2)-10 center(2)+10; -35 -10])
    vf.setRangePosition([-30 10; center(2)-15 center(2)+15; -35 -5])
    
    % Solve for objective origin in restricted volume.
    [origin1, err1] = fminsearch(@(o) objective_origin(o, vf), center);
    origin_obj_res(:, i) = origin1;
    % Solve for objective origin globally.
    vf.setRange(0)
    [origin0, err0] = fminsearch(@(o) objective_origin(o, vf), center);
    origin_obj_res(:, i) = origin1;
    origin_obj_glo(:, i) = origin0;
    
    vfs{i} = vf;
    
    tester3
    
    % Optional visualization.
    vf.plotVector(vf.U_e, 0, '$\vec{u}$');
    hold on
    scatter3(center(1), center(2), center(3), 60, 1, 'filled', ...
        'MarkerFaceColor', 'y')
    hold on
    scatter3(origin1(1), origin1(2), origin1(3), 60, 1, 'filled', ...
        'MarkerFaceColor', 'r')
    hold on
    scatter3(origin0(1), origin0(2), origin0(3), 60, 1, 'filled', ...
        'MarkerFaceColor', 'k')
    legend('$\vec{u}$', 'center', 'restricted objective origin', ...
        'global objective origin')
    pause
    close
end

%%%%%%%%%%%%%%%% progression of objective origin %%%%%%%%%%%%%%%% 
% Restricted.
figure;
scatter3(origin_obj_res(1,:), origin_obj_res(2,:), origin_obj_res(3,:), dotSize, ...
    1: num_frame, 'filled')
% Color corresponds to frame number.
colorbar
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('Objective origin in vortical region over time')

% Global.
figure;
scatter3(origin_obj_glo(1,:), origin_obj_glo(2,:), origin_obj_glo(3,:), dotSize, ...
    1: num_frame, 'filled')
% Color corresponds to frame number.
colorbar
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('Objective origin in global volume over time')

% Compute error associated with the objective origin.
err = zeros(3, num_frame);

% Origin whose impulse value is taken as correct. Currently center of
% measurement volume.
origin_ref = [-15 0 -20]';

% Impulse computed with different origins.
I_ref = zeros(3, num_frame);
% Restricted objective origins.
I_obj_res = zeros(3, num_frame);
I_obj_glo = zeros(3, num_frame);
% Using the center of the feature.
I_cen = zeros(3, num_frame);

for i = 1: num_frame
    vf = vfs{i};
    I_ref(:,i) = vf.impulse(origin_ref, 0);
    I_obj_res(:,i) = vf.impulse(origin_obj_res(:,i), 0);
    I_obj_glo(:,i) = vf.impulse(origin_obj_glo(:,i), 0);
    I_cen(:,i) = vf.impulse([-8 10-1.5*(i-1) -21]', 0);
end

% Consistency of impulse computations. Reference origin used as for
% normalization.
sd_ref = std(I_ref, [], 2);
sdr_obj_res = std(I_obj_res, [], 2) ./ sd_ref;
sdr_obj_glo = std(I_obj_glo, [], 2) ./ sd_ref;
sdr_cen = std(I_cen, [], 2) ./ sd_ref;

mag_err = sqrt(sum(err.^2, 1));

% Plot of error.
figure;
scatter3(origin_obj_res(1,:), origin_obj_res(2,:), origin_obj_res(3,:), dotSize, ...
    mag_err, 'filled')

colorbar
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('Error of objective origin over time')
