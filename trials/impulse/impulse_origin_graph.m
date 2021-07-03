% Prsent the effect of origin selection on the error of impulse
% computation, as shown by a scatter plot of error on various origin
% locations.
%
% Derek Li, July 2021


% Paremeters held constant.
sp = 0.1;
fr = 1;
u0 = 1;

[x, y, z, u, v, w, ~] = hill_vortex_3D(sp, fr, u0, 1);

vf = VelocityField.import_grid_separate(x,y,z,u,v,w);

% Consider stochatic effect of noise introduction.
num_ite = 5;
% Proportional noise.
props = 0: 0.1: 3;

%%%%%%%%%%%%%%%% Graph of error over grid of origins %%%%%%%%%%%%%%%%%%

% Sample origins uniformly from the grid. This is done by making another
% velocity field, which is only used for plotting. The positions of this
% velocity field are the origins.
osp = 1 * ones(1, 3);
oends = [-4 4; -4 4; -4 4];

clear X
[X(:,:,:,1), X(:,:,:,2), X(:,:,:,3)] = meshgrid(oends(1,1): osp(1): oends(1,2), ...
    oends(2,1): osp(2): oends(2,2), oends(3,1): osp(3): oends(3,2));
vfp = VelocityField(X, zeros(size(X)));


% Containers for data across all runs.
% Errors here are mean absolute errors.
err = zeros([vfp.dims 3 num_ite]);
err_box = zeros([vfp.dims 3 num_ite]);
err_gss = zeros([vfp.dims 3 num_ite]);

% Baseline error due to origin selection (and imperfect resolution) without
% noise.
bias_ori = zeros([vfp.dims 3]);

bias_box = zeros([vfp.dims 3]);
bias_gss = zeros([vfp.dims 3]);

% Iterate through grid and compute error profile.
for j = 1: vfp.dims(1)
    for i = 1: vfp.dims(2)
        for k = 1: vfp.dims(3)
            for n = 1: num_ite
                % Run script for impulse error sampling.
                [err(j,i,k,:,n), ~, err_box(j,i,k,:,n), ~, err_gss(j,i,k,:,n), ~, ...
                    bias_box(j,i,k,:), bias_gss(j,i,k,:), bias_ori(j,i,k,:)] = ...
                    impulse_err_stats(vf, props, vfp.X(j,i,k,:), fr, u0);
            end
        end
    end
end

% Average over trials.
err = mean(err, 5);
err_box = mean(err_box, 5);
err_gss = mean(err_gss, 5);

% Take magnitude.
mag_err = squeeze(sqrt(sum(err.^2, 1)));
mag_err_box = squeeze(sqrt(sum(err_box.^2, 1)));
mag_err_gss = squeeze(sqrt(sum(err_gss.^2, 1)));

mag_bias_ori = squeeze(sqrt(sum(bias_ori.^2, 1)));
mag_bias_box = squeeze(sqrt(sum(bias_box.^2, 1)));
mag_bias_gss = squeeze(sqrt(sum(bias_gss.^2, 1)));

%%%%%%%%%%%%%%%% Dimensional Plots %%%%%%%%%%%%%%%%%

% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [2 1 3];
dim_str = {'x', 'y', 'z'};

% Save plots.
img_fdr = strcat('C:\Users\derek\flow\trials\impulse\origin\global\sp=', ...
    string(osp(1)), '\', 'dim=', string(oends(1,2)-oends(1,1)), '\');
mkdir(img_fdr);

for dim = dims
    % Error due to origin selection.
    vfp.plotScalar(squeeze(bias_ori(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ error of resolution and origin'));
    saveas(gcf, strcat(img_fdr, 'bias-ori-', string(dim), '.fig'))
    % Smoother biases.
    vfp.plotScalar(squeeze(bias_box(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ box smoother bias'));
    saveas(gcf, strcat(img_fdr, 'bias-box-', string(dim), '.fig'))
    vfp.plotScalar(squeeze(bias_gss(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ Gaussian smoother bias'));
    saveas(gcf, strcat(img_fdr, 'bias-gss-', string(dim), '.fig'))
    % Noise added.
    vfp.plotScalar(squeeze(err(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ unfiltered error'));
    saveas(gcf, strcat(img_fdr, 'err-unf-', string(dim), '.fig'))
    % Smoothed.
    vfp.plotScalar(squeeze(err_box(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ mean error after box smoothing'));
    saveas(gcf, strcat(img_fdr, 'err-box-', string(dim), '.fig'))
    vfp.plotScalar(squeeze(err_gss(:,:,:,dim)), 0, ...
        strcat('$', dim_str{dim}, '$ mean error after Gaussian smoothing'));
    saveas(gcf, strcat(img_fdr, 'err-gss-', string(dim), '.fig'))
end

% Vector plots of error.
% Error due to origin selection.
vfp.plotVector(bias_ori, 0, ...
    'Error of resolution and origin');
saveas(gcf, strcat(img_fdr, 'bias-ori-v.fig'))
% Smoother biases.
vfp.plotVector(bias_box, 0, ...
    strcat('Box smoother bias'));
saveas(gcf, strcat(img_fdr, 'bias-box-v.fig'))
vfp.plotVector(bias_gss, 0, ...
    strcat('Gaussian smoother bias'));
saveas(gcf, strcat(img_fdr, 'bias-gss-v.fig'))
% Noise added.
vfp.plotVector(err, 0, ...
    strcat('Unfiltered error'));
saveas(gcf, strcat(img_fdr, 'err-unf-v.fig'))
% Smoothed.
vfp.plotVector(err_box, 0, ...
    strcat('Mean error after box smoothing'));
saveas(gcf, strcat(img_fdr, 'err-box-v.fig'))
vfp.plotVector(err_gss, 0, ...
    strcat('Mean error after Gaussian smoothing'));
saveas(gcf, strcat(img_fdr, 'err-gss-v.fig'))


%%%%%%%%%%%%% Magnitude Plots %%%%%%%%%%%%%
% Error due to origin selection.
vfp.plotScalar(mag_bias_ori, 0, ...
    'Magnitude of error of resolution and origin');
saveas(gcf, strcat(img_fdr, 'bias-ori.fig'))
% Smoother biases.
vfp.plotScalar(mag_bias_box, 0, ...
    strcat('Magnitude of box smoother bias'));
saveas(gcf, strcat(img_fdr, 'bias-box.fig'))
vfp.plotScalar(mag_bias_gss, 0, ...
    strcat('Magnitude of Gaussian smoother bias'));
saveas(gcf, strcat(img_fdr, 'bias-gss.fig'))
% Noise added.
vfp.plotScalar(mag_err, 0, ...
    strcat('Magnitude of unfiltered error'));
saveas(gcf, strcat(img_fdr, 'err-unf.fig'))
% Smoothed.
vfp.plotScalar(mag_err_box, 0, ...
    strcat('Magnitude of mean error after box smoothing'));
saveas(gcf, strcat(img_fdr, 'err-box.fig'))
vfp.plotScalar(mag_err_gss, 0, ...
    strcat('Magnitude of mean error after Gaussian smoothing'));
saveas(gcf, strcat(img_fdr, 'err-gss.fig'))
