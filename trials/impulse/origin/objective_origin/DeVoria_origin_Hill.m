% Optimize origin selection by minimizing deviation from the integral objective
% origin constraints computed in objective_origin.m, which is suggested in
% (De Voria 2014).

% Synethetic data set.
spr = 0.05;
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
vf.setRangePosition(range);


% % Use experimental vortex.
% load(sprintf('%s%s', folder, '\data\turbulent_vortex_post.mat'))
% vf = VelocityField.importCmps(x, y, z, u, v, w);% Vortical region.
% vf.setRangePosition([-20 0; -5 25; -35 -5])

% Consider stochatic effect of noise introduction.
num_ite = 20;
% Proportional noise.
props = 0: 0.5: 3;
props_count = size(props, 2);

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

% Record the optimized residuals of the two identities.
res_unf = zeros(props_count, num_ite, 3);
res_box = zeros(props_count, num_ite, 3);
res_gss = zeros(props_count, num_ite, 3);

% Record the optimized residuals of the two identities.
res0_unf = zeros(props_count, num_ite, 3);
res0_box = zeros(props_count, num_ite, 3);
res0_gss = zeros(props_count, num_ite, 3);


% Minimization options.
% min_opt = optimoptions(@fminunc,'Algorithm','trust-region','SpecifyObjectiveGradient',true);
min_opt = optimoptions(@fminunc,'Algorithm','quasi-newton','SpecifyObjectiveGradient',false,'Display','iter');
origin0 = [0 0 0]';
% Initial residual.
res0 = objective_origin_obj(origin0, vf);
[bestorigin, err0, exitflag,~] = fminunc(@(o) objective_origin_obj(o, vf), origin0, min_opt);

% Identify objective origin under different proportions of noise.
for i = 1: props_count
    for j = 1: num_ite
        vf.clearNoise()
        % Add noise.
        N = vf.noise_uniform(props(i)*u0);
        % Identify objective origin without smoothing.
%         % Randomize initial origin guess.
%         origin0 = -2 + 4*rand(3, 1);
%         [origin, res_unf(i,j,2)] = fminunc(@(o) momentum_origin_obj(o, vf, 0), origin0, min_opt);
        [origin, res_unf(i,j,2)] = patternsearch(@(o) momentum_origin_obj(o, vf, 0), origin0);
        % Compute corresponding error.
        err_unf(i,j,:) = vf.impulse(origin, 1) - I0;
        err0_unf(i,j,:) = vf.impulse(origin_ref, 1) - I0;
        origin_unf(i,j,:) = origin;
        % Record residuals.
        res_unf(i,j,1) = momentum_origin_obj(origin, vf, 0);
        res0_unf(i,j,1) = momentum_origin_obj(origin0, vf, 0);
        res0_unf(i,j,2) = objective_origin_obj(origin0, vf);
        
        % Box smoothing and origin identification.
        vf.smoothNoise('box');
%         [origin, res_box(i,j,2)] = fminunc(@(o) momentum_origin_obj(o, vf, 0), origin0, min_opt);
        [origin, res_box(i,j,2)] = patternsearch(@(o) momentum_origin_obj(o, vf, 0), origin0);
        err_box(i,j,:) = vf.impulse(origin, 1) - I0;
        err0_box(i,j,:) = vf.impulse(origin_ref, 1) - I0;
        origin_box(i,j,:) = origin;
        % Record residuals.
        res_box(i,j,1) = momentum_origin_obj(origin, vf, 0);
        res0_box(i,j,1) = momentum_origin_obj(origin0, vf, 0);
        res0_box(i,j,2) = objective_origin_obj(origin0, vf);
        
        % Gaussian smoother and identification.
        vf.setNoise(N)
        vf.smoothNoise('gaussian');
%         [origin, res_gss(i,j,2)] = fminunc(@(o) momentum_origin_obj(o, vf, 0), origin0, min_opt);
        [origin, res_gss(i,j,2)] = patternsearch(@(o) momentum_origin_obj(o, vf, 0), origin0);
        err_gss(i,j,:) = vf.impulse(origin, 1) - I0;
        err0_gss(i,j,:) = vf.impulse(origin_ref, 1) - I0;
        origin_gss(i,j,:) = origin;
        % Record residuals.
        res_gss(i,j,1) = momentum_origin_obj(origin, vf, 0);
        res0_gss(i,j,1) = momentum_origin_obj(origin0, vf, 0);
        res0_gss(i,j,2) = objective_origin_obj(origin0, vf);
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

% Mean error magnitudes.
mean_err_unf = mean(mag_err_unf, 2);
errm_unf_sd = std(mag_err_unf, 0, 2);
mean_err_box = mean(mag_err_box, 2);
errm_box_sd = std(mag_err_box, 0, 2);
mean_err_gss = mean(mag_err_gss, 2);
errm_gss_sd = std(mag_err_gss, 0, 2);

mean_err0_unf = mean(mag_err0_unf, 2);
err0m_unf_sd = std(mag_err0_unf, 0, 2);
mean_err0_box = mean(mag_err0_box, 2);
err0m_box_sd = std(mag_err0_box, 0, 2);
mean_err0_gss = mean(mag_err0_gss, 2);
err0m_gss_sd = std(mag_err0_gss, 0, 2);

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
origins = squeeze(origin_unf(4,:,:));
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
origins = squeeze(origin_box(4,:,:));
scatter3(origins(:,1), origins(:,2), origins(:,3), dot_size, mag_err_box(4,:), 'filled')

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
origins = squeeze(origin_gss(4,:,:));
scatter3(origins(:,1), origins(:,2), origins(:,3), dot_size, mag_err_gss(4,:), 'filled')

cb = colorbar;
cb.Label.String = 'error per $I$';
cb.Label.Interpreter = 'latex';
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('Error of objective origins after Gaussian-filtering', 'FontSize', fsize)

disp('Average error at objective origin:')
fprintf('Unfilteed: %f\n', mean(mag_err_unf(4,:)))
fprintf('Box: %f\n', mean(mag_err_box(4,:)))
fprintf('Gaussian: %f\n', mean(mag_err_gss(4,:)))

% Omitting resolution error and smoothing biases, which are negligible.
disp('Relative errors with respect to reference origin')
fprintf('Unfilteed: %f\n', err_rel_unf(4))
fprintf('Box: %f\n', err_rel_box(4))
fprintf('Gaussian: %f\n', err_rel_gss(4))