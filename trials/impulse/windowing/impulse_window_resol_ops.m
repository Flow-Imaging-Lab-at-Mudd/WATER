function [fres_min_bias, fres_min_err] = ...
    impulse_window_resol_ops(fr, winsize, overlap, fres, origin, err_level, props)
% Identical to impulse_window_resol.m with the allowance of the
% superposition of 0.25, 0.5, 0.75 overlap ratios in one graphs, which are
% also customized.
%
% Derek Li, July 2021.

% Freestream velocity.
u0 = 1;

% Generate range of downsampled spacings for evenly spaced feature
% resolutions.
fres = flip(fres);

sps_count = size(fres, 2);

% Multiple overlaps allowed to be superimposed in graph.
ops_count = length(overlap);

% Spacings for different overlaps.
sps = zeros(sps_count, ops_count);

% Compute initial resolutions to produce the desired given feature
% resolutions after downsampling.
for o = 1: ops_count
    op = round(winsize .* overlap(o));
    op = min([op; winsize-1], [], 1);
    sps(:, o) = 2*fr ./ (2*fres * (winsize-op) + winsize - 1);
end

% Consider stochastic effect.
num_ite = 5;

% Downsampling bias.
dId = zeros(3, sps_count, ops_count);
bias_box = zeros(3, sps_count, ops_count);
bias_gss = zeros(3, sps_count, ops_count);

% Containers for data across all runs.
% Errors here are mean absolute errors.
err = zeros(3, sps_count, num_ite, ops_count);
err_box = zeros(3, sps_count, num_ite, ops_count);
err_gss = zeros(3, sps_count, num_ite, ops_count);
% Standard deviations are also the average over trials. Use SEM instead?
err_sd = zeros(3, sps_count, num_ite, ops_count);
err_sd_box = zeros(3, sps_count, num_ite, ops_count);
err_sd_gss = zeros(3, sps_count, num_ite, ops_count);

for o = 1: ops_count
    for k = 1: sps_count
        % Construct Hill vortex with specified resolution.
        sp = sps(k, o);
        [x, y, z, u, v, w, ~] = hill_vortex_3D(sp, fr, u0, 1);
        vf = VelocityField.import_grid_separate(x,y,z,u,v,w);

        % Subtract freestream velocity to focus on central feature region.
        vf.addVelocity(-vf.U(1,1,1,:))
        % Focus on vortical region.
        vf.setRangePosition(fr*repmat([-1 1], 3, 1))

        for i = 1: num_ite
            % Run script for impulse error sampling.
            [dId(:,k,o), err(:,k,i,o), err_sd(:,k,i,o), err_box(:,k,i,o), ...
                err_sd_box(:,k,i,o), err_gss(:,k,i,o), err_sd_gss(:,k,i,o), ...
                bias_box(:,k,o), bias_gss(:,k,o), ~] = ...
                impulse_err_window_stats(vf, props, origin, fr, u0, winsize, overlap(o));
        end
    end
end

dId = abs(dId);

% Average over trials.
dI = squeeze(mean(err, 3));
dI_sd = squeeze(mean(err_sd, 3));
dI_box = squeeze(mean(err_box, 3));
dI_sd_box = squeeze(mean(err_sd_box, 3));
dI_gss = squeeze(mean(err_gss, 3));
dI_sd_gss = squeeze(mean(err_sd_gss, 3));

% Magnitude of errors.
mag_dId = sqrt(sum(dId.^2, 1));
mag_bias_box = sqrt(sum(bias_box.^2, 1));
mag_bias_gss = sqrt(sum(bias_gss.^2, 1));

mag_dI = sqrt(sum(dI.^2, 1));
mag_dI_sd = sqrt(sum(dI_sd.^2, 1));
mag_dI_box = sqrt(sum(dI_box.^2, 1));
mag_dI_sd_box = sqrt(sum(dI_sd_box.^2, 1));
mag_dI_gss = sqrt(sum(dI_gss.^2, 1));
mag_dI_sd_gss = sqrt(sum(dI_sd_gss.^2, 1));

% Compute the minimum feature resolution required to reduce error beneath
% the desired level.
fres_min_bias = zeros(3, ops_count);
fres_min_err = zeros(3, ops_count);
for o = 1: ops_count
    try
        fres_min_bias(1,o) = fres(end - find(flip(mag_bias_box(:,o)) < err_level, 1) + 1);
    catch
        fres_min_bias(1,o) = -1;
    end
    try
        fres_min_bias(2,o) = fres(end - find(flip(mag_bias_gss(:,o)) < err_level, 1) + 1);
    catch
        fres_min_bias(2,o) = -1;
    end
    try
        fres_min_bias(3,o) = fres(end - find(flip(mag_dId(:,o)) < err_level, 1) + 1);
    catch
        fres_min_bias(3,o) = -1;
    end

    try
        fres_min_err(1,o) = fres(end - find(flip(mag_dI_box(:,o)) < err_level, 1) + 1);
    catch
        fres_min_err(1,o) = -1;
    end
    try
        fres_min_err(2,o) = fres(end - find(flip(mag_dI_gss(:,o)) < err_level, 1) + 1);
    catch
        fres_min_err(2,o) = -1;
    end
    try
        fres_min_err(3,o) = fres(end - find(flip(mag_dI(:,o)) < err_level, 1) + 1);
    catch
        fres_min_err(3,o) = -1;
    end
end


%%%%%%%%%%%%%%%% Dimensional Plots %%%%%%%%%%%%%%%%%
% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [2 1 3];
dim_str = {'x', 'y', 'z'};

% Colors for plotting. From lower to higher overlap, green, blue, red.
smoother_cols = {'g', 'b', 'r'};
cols = repmat([0.7 0.3 0]', 1, 3);


for dim = dims
    % Box smoother bias plot.
    figure;
    % Legend array.
    lgd = cell(1, 2*ops_count);
    
    for o = 1: ops_count
        scatter(fres, dId(dim,:,o), 'ko', 'MarkerFaceColor', cols(o,:), 'LineWidth', 1)
        hold on
        scatter(fres, bias_box(dim,:,o), 'ko', 'MarkerFaceColor', smoother_cols{o}, 'LineWidth', 1)
        hold on
        lgd{2*o-1} = strcat('Downsampling bias, $o=$', {' '}, string(overlap(o)));
        lgd{2*o} = strcat('box smoothed, $o=$', {' '}, string(overlap(o)));
    end
    
    legend(lgd, 'Interpreter', 'latex')
    xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
    ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
    title(strcat('Box $', dim_str{dim}, '$ smoother bias at $r = $', {' '}, string(fr)))
    
    figure;
    for o = 1: ops_count
        scatter(fres, dId(dim,:,o), 'ko', 'MarkerFaceColor', cols(o,:), 'LineWidth', 1)
        hold on
        scatter(fres, bias_gss(dim,:,o), 'ko', 'MarkerFaceColor', smoother_cols{o}, 'LineWidth', 1)
        hold on
        lgd{2*o-1} = strcat('Downsampling bias, $o=$', {' '}, string(overlap(o)));
        lgd{2*o} = strcat('Gaussian smoothed, $o=$', {' '}, string(overlap(o)));
    end

    legend(lgd, 'Interpreter', 'latex')
    xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
    ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
    title(strcat('Gaussian $', dim_str{dim}, '$ smoother bias at $r = $', {' '}, string(fr)))

    % Mean error plot.
    figure;
    for o = 1: ops_count
        errorbar(fres, dI(dim,:,o), dI_sd(dim,:,o), 'ko', ...
            'MarkerFaceColor', cols(o,:), 'LineWidth', 1)
        hold on
        errorbar(fres, dI_box(dim,:,o), dI_sd_box(dim,:,o), 'ko', ...
            'MarkerFaceColor', smoother_cols{o}, 'LineWidth', 1)
        hold on
        
        lgd{2*o-1} = strcat('unfiltered, $o=$', {' '},  string(overlap(o)));
        lgd{2*o} = strcat('box filtered, $o=$', {' '}, string(overlap(o)));
    end
    
    legend(lgd, 'Interpreter', 'latex')
    xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
    ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
    title(strcat('Box $', string(dim_str{dim}), '$ mean error over $\delta u = $', ...
        {' '}, string(props(1)*100), '-', string(props(end)*100), '\% at $r = $', {' '}, string(fr)))
    
    figure;
    for o = 1: ops_count
        errorbar(fres, dI(dim,:,o), dI_sd(dim,:,o), 'ko', ...
            'MarkerFaceColor', cols(o,:), 'LineWidth', 1)
        hold on
        errorbar(fres, dI_gss(dim,:,o), dI_sd_gss(dim,:,o), 'ko', ...
            'MarkerFaceColor', smoother_cols{o}, 'LineWidth', 1)
        hold on
        
        lgd{2*o-1} = strcat('unfiltered, $o=$', {' '}, string(overlap(o)));
        lgd{2*o} = strcat('Gaussian filtered, $o=$', {' '}, string(overlap(o)));
    end

    legend(lgd, 'Interpreter', 'latex')
    xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
    ylabel(strcat('$\left|\frac{\delta I_', string(dim_str{dim}), '}{I}\right|$'))
    title(strcat('Gaussian $', string(dim_str{dim}), '$ mean error over $\delta u = $', ...
        {' '}, string(props(1)*100), '-', string(props(end)*100), '\% at $r = $', {' '}, string(fr)))
end

%%%%%%%%%%%%%%%%%%% Magnitude Plots %%%%%%%%%%%%%%%%%%%%%

% % Smoother bias plot.
% figure;
% % Legend array.
% lgd = cell(1, 3*ops_count);
% 
% for o = 1: ops_count
%     scatter(fres, mag_dId(1,:,o), 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
%     hold on
%     scatter(fres, mag_bias_box(1,:,o), 'ko', 'MarkerFaceColor', 'red', 'LineWidth', 1)
%     hold on
%     scatter(fres, mag_bias_gss(1,:,o), 'ko', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
%     hold on
%     
%     lgd{3*o-2} = strcat('downsampled, $o=$', string(overlap(o)));
%     lgd{3*o-1} = strcat('box smoothed, $o=$', string(overlap(o)));
%     lgd{3*o} = strcat('Gaussian smoothed, $o=$', string(overlap(o)));
% end
% legend(lgd, 'Interpreter', 'latex')
% xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
% ylabel('$\left|\frac{\delta I}{I}\right|$')
% title(strcat('Magnitude of smoother bias at $r = $', {' '}, string(fr)))
% 
% % Mean error plot.
% figure;
% % Legend array.
% lgd = cell(1, 3*ops_count);
% 
% for o = 1: ops_count
%     errorbar(fres, mag_dI(1,:,o), mag_dI_sd(1,:,o), 'ko', 'MarkerFaceColor','black', 'LineWidth', 1)
%     hold on
%     errorbar(fres, mag_dI_box(1,:,o), mag_dI_sd_box(1,:,o), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
%     hold on
%     errorbar(fres, mag_dI_gss(1,:,o), mag_dI_sd_gss(1,:,o), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)
%     hold on
%     
%     lgd{3*o-2} = strcat('unfiltered, $o=$', string(overlap(o)));
%     lgd{3*o-1} = strcat('box filtered, $o=$', string(overlap(o)));
%     lgd{3*o} = strcat('Gaussian filtered, $o=$', string(overlap(o)));
% end
% 
% legend(lgd, 'Interpreter', 'latex')
% xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
% ylabel('$\left|\frac{\delta I}{I}\right|$')
% title(strcat('Mean error magnitude over $\delta u = $', ...
%         string(props(1)*100), '-', string(props(end)*100), '\% at $r = $', {' '}, string(fr)))


