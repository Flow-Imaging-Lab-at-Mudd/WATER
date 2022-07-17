% Prsent the effect of origin selection on the error of impulse
% computation, as shown by a scatter plot of error on various origin
% locations.
%
% Currently, the experimental data and synthetic data are processed
% differently, which should be more effectively organized.

% Whether figures generated are to be automatically saved.
savefig = false;

% Paremeters of vortex.
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
vf.setRangePosition(fr*repmat([-1 1], 3, 1))

% % Experimental data set.
% load(sprintf('%s%s', rootFolder, '\data\turbulent_vortex_post.mat'))
% vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
% % % Vortical region.
% % vf.setRangePosition([-20 0; -5 25; -35 -5])
% % Central location.
% center = [-8 10 -21]';

% Consider stochatic effect of noise introduction.
num_ite = 5;
% Proportional noise. For biases to be properly computed, must include 0.
props = [0 1.5];
props_count = size(props, 2);


% Theoretical momentum.
I0 = Hill_Impulse(vf.fluid.density, vf.scale.len, fr, u0);
i0 = I0(3);

% % Original momentum computed from central location.
% I0 = vf.impulse(center, 0);
% i0 = norm(I0);

%%%%%%%%%%%%%%%% Graph of error over grid of origins %%%%%%%%%%%%%%%%%%

% Sample origins uniformly from the grid. This is done by making another
% velocity field, which is only used for plotting. The positions of this
% velocity field are the origins.
osp = 0.25 * ones(1, 3);
oends = 1.25*repmat([-1 1], 3, 1);

% % Experimental.
% oends = [-30 20; -15 35; -45 -5];
% osp = 10 * ones(1, 3);

clear X
[X(:,:,:,1), X(:,:,:,2), X(:,:,:,3)] = meshgrid(oends(1,1): osp(1): oends(1,2), ...
    oends(2,1): osp(2): oends(2,2), oends(3,1): osp(3): oends(3,2));
vfp = VelocityField(X, zeros(size(X)));


% Containers for data across all runs.
% Errors here are mean absolute errors.
err = zeros([vfp.dims 3 props_count num_ite]);
err_box = zeros([vfp.dims 3 props_count num_ite]);
err_gss = zeros([vfp.dims 3 props_count num_ite]);

% Iterate through grid and compute error profile.
% Plot energy estimation error for small and large values of noise.
for n = 1: num_ite
    % This section parallels impulse_err_run.m
    for p = 1: props_count
        vf.clearNoise()
        N = vf.noise_uniform(props(p)*u0);
        for j = 1: vfp.dims(1)
            for i = 1: vfp.dims(2)
                for k = 1: vfp.dims(3)
                    err(j,i,k,:,p,n) = vf.impulse(squeeze(vfp.X(j,i,k,:)), 1) - I0;
                end
            end
        end
        % Result with box smoothing.
        vf.smoothNoise('box');
        for j = 1: vfp.dims(1)
            for i = 1: vfp.dims(2)
                for k = 1: vfp.dims(3)
                    err_box(j,i,k,:,p,n) = vf.impulse(squeeze(vfp.X(j,i,k,:)), 1) - I0;
                end
            end
        end
        % Reset and smooth with gaussian filter.
        vf.setNoise(N)
        vf.smoothNoise('gaussian');
        for j = 1: vfp.dims(1)
            for i = 1: vfp.dims(2)
                for k = 1: vfp.dims(3)
                    err_gss(j,i,k,:,p,n) = vf.impulse(squeeze(vfp.X(j,i,k,:)), 1) - I0;
                end
            end
        end
    end
end

% Normalize absolute error.
err0 = err / i0;
err0_box = err_box / i0;
err0_gss = err_gss / i0;

% Extract bias from a 0 noise level.
bias_ori = squeeze(err(:,:,:,:,1,1));
bias_box = squeeze(err_box(:,:,:,:,1,1));
bias_gss = squeeze(err_gss(:,:,:,:,1,1));

% Compute bias magnitude.
mag_bias_ori = squeeze(sqrt(sum(bias_ori.^2, 4)));
mag_bias_box = squeeze(sqrt(sum(bias_box.^2, 4)));
mag_bias_gss = squeeze(sqrt(sum(bias_gss.^2, 4)));

% Mean error magnitudes.
mag_err = squeeze(mean(sqrt(sum(err.^2, 4)), 6));
mag_err_box = squeeze(mean(sqrt(sum(err_box.^2, 4)), 6));
mag_err_gss = squeeze(mean(sqrt(sum(err_gss.^2, 4)), 6));

% Average over trials.
err = mean(abs(err0), 6);
err_box = mean(abs(err0_box), 6);
err_gss = mean(abs(err0_gss), 6);

%%%%%%%%%%%%%%%% Dimensional Plots %%%%%%%%%%%%%%%%%

% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [];
dim_str = {'x', 'y', 'z'};
% Font size for printing titles and axis labels.
vfp.setFontSize(11);

% Save plots.
savePlot = 0;
if savePlot
    img_fdr = strcat(rootFolder, '\trials\impulse\origin\global\du=', ...
        props(2)*100, '\%', '\');
    if ~isfolder(img_fdr)
        mkdir(img_fdr);
    end
end

for dim = dims
    % Error due to origin selection.
    vfp.plotScalar(squeeze(bias_ori(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ error of resolution and origin'));
    if savePlot
        saveas(gcf, strcat(img_fdr, 'bias-ori-', string(dim), '.fig'))
    end
    % Smoother biases.
    vfp.plotScalar(squeeze(bias_box(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ box smoother bias'));
    if savePlot
        saveas(gcf, strcat(img_fdr, 'bias-box-', string(dim), '.fig'))
    end
    vfp.plotScalar(squeeze(bias_gss(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ Gaussian smoother bias'));
    if savePlot
        saveas(gcf, strcat(img_fdr, 'bias-gss-', string(dim), '.fig'))
    end
    % Noise added.
    vfp.plotScalar(squeeze(err(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ unfiltered error'));
    if savePlot
        saveas(gcf, strcat(img_fdr, 'err-unf-', string(dim), '.fig'))
    end
    % Smoothed.
    vfp.plotScalar(squeeze(err_box(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ mean error after box smoothing'));
    if savePlot
        saveas(gcf, strcat(img_fdr, 'err-box-', string(dim), '.fig'))
    end
    vfp.plotScalar(squeeze(err_gss(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ mean error after Gaussian smoothing'));
    if savePlot
        saveas(gcf, strcat(img_fdr, 'err-gss-', string(dim), '.fig'))
    end
end

% % Vector plots of error.
% % Error due to origin selection.
% vfp.plotVector(bias_ori, 0, ...
%     'Error of resolution and origin');
% saveas(gcf, strcat(img_fdr, 'bias-ori-v.fig'))
% % Smoother biases.
% vfp.plotVector(bias_box, 0, ...
%     strcat('Box smoother bias'));
% saveas(gcf, strcat(img_fdr, 'bias-box-v.fig'))
% vfp.plotVector(bias_gss, 0, ...
%     strcat('Gaussian smoother bias'));
% saveas(gcf, strcat(img_fdr, 'bias-gss-v.fig'))
% % Noise added.
% vfp.plotVector(err, 0, ...
%     strcat('Unfiltered error'));
% saveas(gcf, strcat(img_fdr, 'err-unf-v.fig'))
% % Smoothed.
% vfp.plotVector(err_box, 0, ...
%     strcat('Mean error after box smoothing'));
% saveas(gcf, strcat(img_fdr, 'err-box-v.fig'))
% vfp.plotVector(err_gss, 0, ...
%     strcat('Mean error after Gaussian smoothing'));
% saveas(gcf, strcat(img_fdr, 'err-gss-v.fig'))


%%%%%%%%%%%%% Magnitude Plots %%%%%%%%%%%%%
% Error due to origin selection.
vfp.plotScalar(mag_bias_ori, 0, ...
    'Magnitude of error of resolution and origin');
if savePlot
    saveas(gcf, strcat(img_fdr, 'bias-ori.fig'))
end
% Smoother biases.
vfp.plotScalar(mag_bias_box, 0, ...
    strcat('Magnitude of box smoother bias'));
if savePlot
    saveas(gcf, strcat(img_fdr, 'bias-box.fig'))
end
vfp.plotScalar(mag_bias_gss, 0, ...
    strcat('Magnitude of Gaussian smoother bias'));
if savePlot
    saveas(gcf, strcat(img_fdr, 'bias-gss.fig'))
end

% Noise added.
vfp.plotScalar(mag_err(:,:,:,2), 0, ...
    strcat('Magnitude of unfiltered error'));
if savePlot
    saveas(gcf, strcat(img_fdr, 'err-unf.fig'))
end
% Smoothed.
vfp.plotScalar(mag_err_box(:,:,:,2), 0, ...
    strcat('Magnitude of mean error after box smoothing'));
if savePlot
    saveas(gcf, strcat(img_fdr, 'err-box.fig'))
end
vfp.plotScalar(mag_err_gss(:,:,:,2), 0, ...
    strcat('Magnitude of mean error after Gaussian smoothing'));
if savePlot
    saveas(gcf, strcat(img_fdr, 'err-gss.fig'))
end
% 
% % % Save workspace.
% % save(strcat(img_fdr, 'data', string(size(X, 2)), '.mat'))
% 
disp('Minimum error on grid')
disp(min(mag_err(:)))

disp('Mean error sampled from grid:')
fprintf('Unfilteed: %f\n', mean(mag_err, 'all'))
fprintf('Box: %f\n', mean(mag_err_box, 'all'))
fprintf('Gaussian: %f\n', mean(mag_err_gss, 'all'))

% % Optimal origin throgh minimization on THE LAST error profile. Used when a
% % single error profile is considered.
% [origin_min, err_min] = optimal_origin(vf, I0, [0 0 0]', oends);
% disp('Minimum Error within Grid')
% disp(err_min)
% origin_min
