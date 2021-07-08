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

% Consider stochatic effect of noise introduction.
num_ite = 5;
% Proportional noise.
props = 0: 0.5: 3;
props_count = size(props, 2);

% Set constant maximal magnitude of noise.
u_mean = vf.meanSpeed(0, 0);

% Theoretical momentum.
I0 = vf.fluid.density*[0 2*pi*fr^3*u0*vf.scale.len^4 0]';
i0 = I0(2);

% Objective origins identified per noise level per trial. For no filtering,
% box filtering, and Gaussian filtering.
origin_unf = zeros(props_count, num_ite, 3);
err_unf = zeros(props_count, num_ite, 3);
origin_box = zeros(props_count, num_ite, 3);
err_box = zeros(props_count, num_ite, 3);
origin_gss = zeros(props_count, num_ite, 3);
err_gss = zeros(props_count, num_ite, 3);

% Identify objective origin under different proportions of noise.
for i = 1: props_count
    for j = 1: num_ite
        vf.clearNoise();
        % Add noise.
        N = vf.noise_uniform(props(i)*u_mean);
        % Identify objective origin without smoothing.
        origin = fminsearch(@(o) objective_origin(o, vf), [-10 10 15]');
        % Compute corresponding error.
        err_unf(i,j,:) = vf.impulse(1, origin) - I0;
        origin_unf(i,j,:) = origin;
        % Box smoothing and origin identification.
        vf.smoothNoise('box');
        origin = fminsearch(@(o) objective_origin(o, vf), [-10 10 15]');
        err_box(i,j,:) = vf.impulse(1, origin) - I0;
        origin_box(i,j,:) = origin;
        % Gaussian smoother and identification.
        vf.setNoise(N)
        vf.smootheNoise('gaussian');
        origin = fminsearch(@(o) objective_origin(o, vf), [-10 10 15]');
        err_gss(i,j,:) = vf.impulse(1, origin) - I0;
        origin_box(i,j,:) = origin;
    end
end

err_unf = err_unf / i0;
mag_err_unf = sqrt(sum(err_unf.^2, 3));

err_box = err_box / i0;
mag_err_box = sqrt(sum(err_box.^2, 3));

err_gss = err_gss / i0;
mag_err_gss = sqrt(sum(err_gss.^2, 3));

%%%%%%%%%%%%%%%%%% Scatter plot of objective origins %%%%%%%%%%%%%%%%%%%
% Unfiltered.
origins = reshape(origin_unf, [], 3);
err = mag_err_unf(:);

    
    
    
