function [dI, dI_box, dI_gss, dI0, bias_box, bias_gss, ...
    di, di_box, di_gss, di0, mag_bias_box, mag_bias_gss, ...
    dI_sd, dI_sd_box, dI_sd_gss, di_sd, di_sd_box, di_sd_gss, axes] = ...
    impulse_winsize(vf, I0, origin, props, winsizes, overlap, num_ite, display_plots, Ih)
% Vary the window size used in downsampling and present its effect on error.
%
% April, 2022

font = 'Arial';
fontSize = 8;

if ~isvector(winsizes)
    error('Uniform windows expected!')
end

if ~exist('Ih', 'var')
    Ih = @(vf, with_noise) vf.impulse(origin, with_noise);
end
    
% Assume uniform windows in x, y, z.
win_count = length(winsizes);

% Containers for error data at different window sizes.
dI = zeros(3, win_count);
dI_sd = zeros(3, win_count);
dI0 = zeros(3, win_count);
dI_box = zeros(3, win_count);
dI_sd_box = zeros(3, win_count);
dI_gss = zeros(3, win_count);
dI_sd_gss = zeros(3, win_count);
bias_box = zeros(3, win_count);
bias_gss = zeros(3, win_count);

% Magnitude errors.
di = zeros(1, win_count);
di_sd = zeros(1, win_count);
di0 = zeros(1, win_count);
di_box = zeros(1, win_count);
di_sd_box = zeros(1, win_count);
mag_bias_box = zeros(1, win_count);
di_gss = zeros(1, win_count);
di_sd_gss = zeros(1, win_count);
mag_bias_gss = zeros(1, win_count);

for k = 1: win_count
    disp(['Window size ' num2str(k)])
    [dI(:,k), dI_box(:,k), dI_gss(:,k), dI0(:,k), bias_box(:,k), bias_gss(:,k), ...
        dI_sd(:,k), dI_sd_box(:,k), dI_sd_gss(:,k), ...
        di(:,k), di_box(:,k), di_gss(:,k), di0(:,k), mag_bias_box(:,k), mag_bias_gss(:,k), ...
        di_sd(:,k), di_sd_box(:,k), di_sd_gss(:,k), ~] = ...
        impulse_err_run_constN(vf, props, origin, I0, num_ite, [winsizes(k), overlap], {}, Ih);
end

%%%%%%%%%%% Dimensional Plots of Error %%%%%%%%%%%%%
if ~exist('display_plots', 'var') || ((isinteger(display_plots)||islogical(display_plots)) && ~display_plots)
    return
end
axes = {};
% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [3];
dim_str = {'x', 'y', 'z'};

if ismember('dim', display_plots)
    for dim = dims
        % Plot of baseline resolution error.
        if ismember('win', display_plots)
%             figure;
            axes{end+1} = scatter(winsizes, dI0(dim,:), 'k', 'filled');
            xticks(winsizes)
            xlabel('Window size (voxels)','fontName',font,'fontSize',fontSize,'interpreter','none')
            ylabel('$\frac{\delta I}{|\vec{I}|}$', 'HorizontalAlignment', 'right', 'Rotation', 0,'fontSize',1.5*fontSize)
            title('Impulse windowing resolution error in $\hat{z}$')
            % Log plot.
            ax = gca;
            %ax.XScale = 'log';
            xlim([winsizes(1)/2 winsizes(end)*2])
        end
        
        % Plot of windowing error filter errors.
        if ismember('resol', display_plots)
%             figure;
            axes{end+1} = scatter(winsizes, dI0(dim,:), 'k', 'filled');
%             hold on
%             scatter(winsizes, bias_box(dim,:), 'r', 'filled', 'Marker', 's')
            hold on
            scatter(winsizes, bias_gss(dim,:), 'b', 'filled', 'Marker', '^')
            
            xticks(winsizes)
            legend({'Unfiltered', 'Gaussian'},'fontName',font,'fontSize',fontSize,'interpreter','none','location','southwest')
%             legend({'unfiltered', 'box', 'Gaussian'})
            xlabel('Window size (voxels)','fontName',font,'fontSize',fontSize,'interpreter','none')
            ylabel('$\frac{\delta I}{|\vec{I}|}$', 'HorizontalAlignment', 'right', 'Rotation', 0,'fontSize',1.5*fontSize)
            title('Impulse windowing resolution error in $\hat{z}$')
            % Log plot.
            ax = gca;
            %ax.XScale = 'log';
            xlim([winsizes(1)/2 winsizes(end)*2])
        end
        
        % Plot of mean errors.
        if ismember('mean', display_plots)
            %     figure;
            axes{end+1} = errorbar(winsizes, dI(dim,:), dI_sd(dim,:), 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1);
            hold on
            errorbar(winsizes, dI_box(dim,:), dI_sd_box(dim,:), 'Marker', '+', 'MarkerEdgeColor', 'red', 'LineStyle', 'none')
            hold on
            errorbar(winsizes, dI_gss(dim,:), dI_sd_gss(dim,:), 'Marker', 'x', 'MarkerEdgeColor', 'blue', 'LineStyle', 'none')
            
            xticks(winsizes)
            legend({'Unfiltered', 'Gaussian'})
%             legend({'unfiltered', 'box', 'Gaussian'})
            xlabel('Window size (voxels)','fontName',font,'fontSize',fontSize,'interpreter','none')
            ylabel('$\frac{\delta I}{|\vec{I}|}$', 'HorizontalAlignment', 'right', 'Rotation', 0,'fontSize',1.5*fontSize)
            title('Impulse error in $\hat{z}$ under noise')
            % Log plot.
            ax = gca;
            %ax.XScale = 'log';
        end
    end
end

%%%%%%%%%%% Magnitude plots %%%%%%%%%%%%%%
% % Plot of windowing error.
% figure;
% scatter(winsizes, di0, 'k', 'filled')
% xticks(winsizes)
% xlabel('Window size')
% ylabel('Normalized error')
% title('Windowing resolution error of impulse')
% 
% % Plot of bias.
% figure;
% scatter(winsizes, di0, 'k', 'filled')
% hold on
% scatter(winsizes, mag_bias_box, 'r', 'filled')
% hold on
% scatter(winsizes, mag_bias_gss, 'b', 'filled')

% xticks(winsizes)
% legend({'windowing', 'box', 'Gaussian'})
% xlabel('Window size')
% ylabel('Normalized error')
% title('Baseline impulse resolution error of filters')

% Plot of mean absolute errors.
if ismember('mag', display_plots)
%     figure;
    errorbar(winsizes, di, di_sd, 'ko', 'MarkerFaceColor', 'black', 'LineWidth', 1)
%     hold on
%     errorbar(winsizes, di_box, di_sd_box, 'Marker', 's', 'MarkerFaceColor', 'red', 'Color', 'red', 'LineWidth', 1, 'LineStyle', 'none')
    hold on
    axes{end+1} = errorbar(winsizes, di_gss, di_sd_gss, 'Marker', '^', 'MarkerFaceColor', 'blue', 'Color', 'blue', 'LineWidth', 1, 'LineStyle', 'none');
    xticks(winsizes)
    legend({'Unfiltered', 'Gaussian'})
%     legend({'unfiltered', 'box', 'Gaussian'})
    xlabel('Window size')
    ylabel('$\frac{|\delta \vec{I}|}{|\vec{I}|}$', 'HorizontalAlignment', 'right', 'Rotation', 0)
    title('Impulse error magnitude under noise')
    % Log plot.
    ax = gca;
    %ax.XScale = 'log';
    xlim([winsizes(1)/2 winsizes(end)*2])
end
end
