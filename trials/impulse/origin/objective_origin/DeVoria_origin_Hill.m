% Optimize origin selection by minimizing deviation from the integral objective
% origin constraints computed in objective_origin.m, which is suggested in
% (De Voria 2014).

% Synethetic data set.
spr = 0.1;
l = 1;
vr = 1;
% Feature radius.
fr = l*vr;
u0 = 1;
% Remove free stream.
[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, u0, 1);

vf = VelocityField.importCmps(x, y, z, u, v, w);
% Zoom in on vortical region.
range = fr*repmat([-1 1], 3, 1);
range(3,1)=0;
vf.setRangePosition(range);


% % Use experimental vortex.
% load(sprintf('%s%s', folder, '\data\turbulent_vortex_post.mat'))
% vf = VelocityField.importCmps(x, y, z, u, v, w);% Vortical region.
% vf.setRangePosition([-20 0; -5 25; -35 -5])

% Consider stochatic effect of noise introduction.
num_ite = 5;
% Proportional noise.
props = [0 2];
props_count = size(props, 2);

% Set constant maximal magnitude of noise.
u_mean = vf.meanSpeed(0, 0);

% Theoretical momentum.
I0 = Hill_Impulse(vf.fluid.density, vf.scale.len, vr, u0, vr);
i0 = I0(3);

% % Momentum computed without noise.
% I0 = vf.impulse([0 0 0]', 0);
% i0 = norm(I0);

% Objective origins identified per noise level per trial. For no filtering,
% box filtering, and Gaussian filtering.
origin_unf = zeros(props_count, num_ite, 3);
err_unf = zeros(props_count, num_ite, 3);
origin_box = zeros(props_count, num_ite, 3);
err_box = zeros(props_count, num_ite, 3);
origin_gss = zeros(props_count, num_ite, 3);
err_gss = zeros(props_count, num_ite, 3);

% Central origin used for comparison.
origin_ref = [0 0 0]';
err0_unf = zeros(props_count, num_ite, 3);
err0_box = zeros(props_count, num_ite, 3);
err0_gss = zeros(props_count, num_ite, 3);

% Minimization options.
%min_opt = optimoptions(@fminunc,'Algorithm','trust-region','SpecifyObjectiveGradient',true);
min_opt = optimoptions(@fminunc,'Algorithm','quasi-newton','SpecifyObjectiveGradient',false,'Display','iter');
origin0 = [0 0 0]';
err0=objective_origin_obj(origin0, vf); 
[bestorigin,err0,exitflag,~] = fminunc(@(o) objective_origin_obj(o, vf), origin0, min_opt);

% Identify objective origin under different proportions of noise.
for i = 1: props_count
    for j = 1: num_ite
        vf.clearNoise()
        % Add noise.
        N = vf.noise_uniform(props(i)*u_mean);
        % Identify objective origin without smoothing.
%         % Randomize initial origin guess.
%         origin0 = -2 + 4*rand(3, 1);
        origin = fminunc(@(o) objective_origin_obj(o, vf), origin0, min_opt);
        % Compute corresponding error.
        err_unf(i,j,:) = vf.impulse(origin, 1) - I0;
        err0_unf(i,j,:) = vf.impulse(origin_ref, 1) - I0;
        origin_unf(i,j,:) = origin;
        % Box smoothing and origin identification.
        vf.smoothNoise('box');
        origin = fminunc(@(o) objective_origin_obj(o, vf), origin0, min_opt);
        err_box(i,j,:) = vf.impulse(origin, 1) - I0;
        err0_box(i,j,:) = vf.impulse(origin_ref, 1) - I0;
        origin_box(i,j,:) = origin;
        % Gaussian smoother and identification.
        vf.setNoise(N)
        vf.smoothNoise('gaussian');
        origin = fminunc(@(o) objective_origin_obj(o, vf), origin0, min_opt);
        err_gss(i,j,:) = vf.impulse(origin, 1) - I0;
        err0_gss(i,j,:) = vf.impulse(origin_ref, 1) - I0;
        origin_gss(i,j,:) = origin;
    end
end

err_unf = err_unf / i0;
mag_err_unf = sqrt(sum(err_unf.^2, 3));
err0_unf = err0_unf / i0;
mag_err0_unf = sqrt(sum(err0_unf.^2, 3));

err_box = err_box / i0;
mag_err_box = sqrt(sum(err_box.^2, 3));
err0_box = err0_box / i0;
mag_err0_box = sqrt(sum(err0_box.^2, 3));

err_gss = err_gss / i0;
mag_err_gss = sqrt(sum(err_gss.^2, 3));
err0_gss = err0_gss / i0;
mag_err0_gss = sqrt(sum(err0_gss.^2, 3));

% Mean error profiles.
mean_err_unf = mean(mag_err_unf, 2);
mean_err_box = mean(mag_err_box, 2);
mean_err_gss = mean(mag_err_gss, 2);

mean_err0_unf = mean(mag_err0_unf, 2);
mean_err0_box = mean(mag_err0_box, 2);
mean_err0_gss = mean(mag_err0_gss, 2);

% Relative errors.
err_rel_unf = mean_err_unf - mean_err0_unf;
err_rel_box = mean_err_box - mean_err0_box;
err_rel_gss = mean_err_gss - mean_err0_gss;


%%%%%%%%%%%%%%%%%% Scatter plot of objective origins %%%%%%%%%%%%%%%%%%%
% Font size for title.
fsize = 15;

% Unfiltered.
figure;
dot_size = 20;
origins = squeeze(origin_unf(2,:,:));
scatter3(origins(:,1), origins(:,2), origins(:,3), dot_size, mag_err_unf(2,:), 'filled')

cb = colorbar;
cb.Label.String = 'error per $I$';
cb.Label.Interpreter = 'latex';
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('Error of objective origins without filtering', 'FontSize', fsize)

% Box filtered.
figure;
dot_size = 20;
origins = squeeze(origin_box(2,:,:));
scatter3(origins(:,1), origins(:,2), origins(:,3), dot_size, mag_err_box(2,:), 'filled')

cb = colorbar;
cb.Label.String = 'error per $I$';
cb.Label.Interpreter = 'latex';
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('Error of objective origins after box-filtering', 'FontSize', fsize)

% Gaussian Filtered.
figure;
dot_size = 20;
origins = squeeze(origin_gss(2,:,:));
scatter3(origins(:,1), origins(:,2), origins(:,3), dot_size, mag_err_gss(2,:), 'filled')

cb = colorbar;
cb.Label.String = 'error per $I$';
cb.Label.Interpreter = 'latex';
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('Error of objective origins after Gaussian-filtering', 'FontSize', fsize)

disp('Average error at objective origin:')
fprintf('Unfilteed: %f\n', mean(mag_err_unf(2,:)))
fprintf('Box: %f\n', mean(mag_err_box(2,:)))
fprintf('Gaussian: %f\n', mean(mag_err_gss(2,:)))

% Omitting resolution error and smoothing biases, which are negligible.
disp('Relative errors with respect to reference origin')
fprintf('Unfilteed: %f\n', err_rel_unf(2))
fprintf('Box: %f\n', err_rel_box(2))
fprintf('Gaussian: %f\n', err_rel_gss(2))