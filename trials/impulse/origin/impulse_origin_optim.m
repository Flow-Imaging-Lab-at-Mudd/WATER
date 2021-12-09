% Identify origin of (locally) minimum error over a number of random error
% trials with potentially varying levels of noise, in company to
% impulse_origin_graph.m
%
% Derek Li, July 2021

% Synthetic data set.
% Paremeters held constant.
sp = 0.1;
fr = 1;
u0 = 1;
[x, y, z, u, v, w, ~] = Hill_Vortex(sp, fr, u0, 1, 1);

vf = VelocityField.importCmps(x, y, z, u, v, w);
% Zoom in on vortical region.
vf.setRangePosition(fr*repmat([-1 1], 3, 1))
% Theoretical impulse.
I0 = HillImpulse(vf.fluid.density, vf.scale.len, fr, u0);

% % Use experimental vortex.
% load(sprintf('%s%s', folder, '\data\turbulent_vortex_post.mat'))
% vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
% % Vortical region.
% vf.setRangePosition([-20 0; -5 25; -35 -5])

% Proportional noise.
props = 0: 0.5: 3;

% Sample origins uniformly from the grid. This is done by making another
% velocity field, which has no velocity content. The positions of this
% velocity field are the origins.
% Synthetic.
osp = 2 * ones(1, 3);
oends = [-2 2; -2 2; -2 2];

% % Experimental.
% oends = [-30 20; -15 35; -45 -5];
% osp = 15 * ones(1, 3);

clear X
[X(:,:,:,1), X(:,:,:,2), X(:,:,:,3)] = meshgrid(oends(1,1): osp(1): oends(1,2), ...
    oends(2,1): osp(2): oends(2,2), oends(3,1): osp(3): oends(3,2));
vfp = VelocityField(X, zeros(size(X)));

% Iterations of error trials to assess average error.
num_ite = 5;

% Container for optimal solutions found upon different initial origin
% guesses.
origin = zeros([vfp.dims 3]);
err = zeros(vfp.dims);

% Boundaries for pattern search minimization algorithm.
ps_range = 2*repmat([-1 1], 3, 1)';
ps_opt = optimoptions(@patternsearch,'Display','iter');

% origin0 = [0 0 0]';
origin0 = [-8 10 -21]';

% Minimize with different initial guesses uniformly distributed over the
% grid.
for j = 1: vfp.dims(1)
    for i = 1: vfp.dims(2)
        for k = 1: vfp.dims(3)
            origin0 = squeeze(vfp.X(j,i,k,:));
            [origin(j,i,k,:), err(j,i,k)] = patternsearch(...
                @(o) min_err_origin(o, 'none', num_ite, vf, props, I0), ...
                origin0, [],[],[],[], ps_range(1,:), ps_range(2,:), ps_opt);
        end
    end
end

figure;
dot_size = 20;
origins = reshape(origin, [], 3);
scatter3(origins(:,1), origins(:,2), origins(:,3), dot_size, err(:), 'filled')

colorbar
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('Error of optimal origin')
