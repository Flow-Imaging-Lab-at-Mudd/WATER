% Visualization of the objective function and its gradient.
%
% Derek Li, July 2021

% Paremeters held constant.
sp = 0.1;
fr = 1;
u0 = 1;

[x, y, z, u, v, w, ~] = hill_vortex_3D(sp, fr, u0, 1);

vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
% Subtract freestream.
vf.addVelocity(-vf.U(1,1,1,:))
% Zoom in on vortical region.
vf.setRangePosition(fr*repmat([-1 1], 3, 1))

% Introduce noise to the system.
% vf.noise_uniform(3*vf.meanSpeed(0, 0));

% Range of points to plot the objective function, uniform in 3 dimensions.
range = -1: 0.2: 1;
count = size(range, 2);

clear X
[X(:,:,:,1), X(:,:,:,2), X(:,:,:,3)] = meshgrid(-2:0.5:2, -2:0.5:2, -2:0.5:2);
vfp = VelocityField(X, zeros(size(X)));

% Error from the objective function.
err = zeros(count, count, count);
grad = zeros(count, count, count, 3);

% Iterate through grid.
for j = 1: vfp.dims(1)
    for i = 1: vfp.dims(2)
        for k = 1: vfp.dims(3)
            % Run script for impulse error sampling.
            [err(j,i,k), grad(j,i,k,:)] = ...
                objective_origin(squeeze(vfp.X(j,i,k,:)), vf);
        end
    end
end

% Plot error and negative gradient.
vfp.plotScalar(err, 0, 'Error from objective origin')
vfp.plotVector(-grad, 0, 'Descent direction')
