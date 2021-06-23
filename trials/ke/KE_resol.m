% Uniform resolutions considered.
resol = 0.01: 0.01: 1;
resol_count = length(resol);

% Constant feature radius used while varying global resolution.
fr = 1;

% Downsampling parameters.
winsize = 8;
overlap = 0.75;

% Randomly sample effective regions.
num_ite = 1;

% Containers for data across all runs.
% Errors here are mean absolute errors.
err_box = zeros(resol_count, num_ite);
err_gss = zeros(resol_count, num_ite);

bias_box = zeros(resol_count, num_ite);
bias_gss = zeros(resol_count, num_ite);

% alpha_box = zeros(resol_count, num_ite);
% alpha_gss = zeros(resol_count, num_ite);

% Downsampling bias.
dKd = zeros(resol_count, num_ite);

% Global spacing of downsampled data.
dsp = min((1-overlap)*winsize, winsize-1)*resol;

% Flag for under-resolution.
under_resolved = false;

% Introduce noise proportionally.
props = 0: 0.1: 1.5;

for k = 1: resol_count
    % Construct Hill vortex with specified resolution.
    sp = resol(k);
    [x, y, z, u, v, w, Mag] = hill_vortex_3D(sp, fr, 1, 1);
    vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
    
    % Subtract freestream velocity to focus on central feature region.
    vf.addVelocity(-vf.U(1,1,1,:))
    
    
    for i = 1: num_ite
        % This section mirrors KE_uniform.
        % Run script for KE error samrpling.
%         [err_box(k, i), err_gss(k, i), bias_box(k), ...
%             bias_gss(k), alpha_box(k, i), alpha_gss(k, i)] = ...
%             KE_mean_err_run(vf, props, fr);
        % Skip if under-resolved.
        try
        % Incorporate windowing.
        [~, dKd(k,i),  ~, ~, bias_box(k,i), bias_gss(k,i), ...
            err_box(k,i), err_gss(k,i)] = ...
            KE_err_window_run(vf, winsize, overlap, props, fr);
        catch
            % Lower resolutions impossible.
            dsp = dsp(1:k-1);
            dKd = dKd(1:k-1, :);
            bias_box = bias_box(1:k-1, :);
            bias_gss = bias_gss(1:k-1, :);
            err_box = err_box(1:k-1, :);
            err_gss = err_gss(1:k-1, :);
            % Set flag for under-resolution to exit iterations.
            under_resolved = true;
            break
        end
    end
    
    if under_resolved
        break;
    end
end

mean_err_box = mean(err_box, 2);
mean_err_gss = mean(err_gss, 2);

smoother_bias_box = mean(bias_box, 2);
smoother_bias_gss = mean(bias_gss, 2);

% amp_box = mean(alpha_box, 2);
% amp_gss = mean(alpha_gss, 2);

% Downsampling bias.
dbias = mean(dKd, 2);

figure;
% errorbar(resol, mean_err_box, std(err_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth',1)
% hold on
% errorbar(resol, mean_err_gss, std(err_gss, 0, 2), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)

scatter(fr./dsp, mean_err_box, 'ko', 'MarkerFaceColor','red', 'LineWidth',1)
hold on
scatter(fr./dsp, mean_err_gss, 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)

% Mean filter error plot.
legend('box-filtered', 'Gaussian-filtered')
xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
ylabel('$\bar{\left|\frac{\delta K}{K}\right|}$')
title(strcat('Mean Error of Smoother over $\delta u = $', {' '}, ...
    string(props(1)), '-', string(props(end)*100), '\% $\bar{u}$', ...
    ' at $r = $', {' '}, string(fr)))

% Folder for plots.
img_fdr = strcat('C:\Users\derek\flow\trials\ke\resol\');
mkdir(strcat(img_fdr, 'mean-error\'));
mkdir(strcat(img_fdr, 'smoother-bias'));
mkdir(strcat(img_fdr, 'amp-coeff'))
% saveas(gcf, strcat(img_fdr, '\mean-error\err-', string(fr), 'r', '.jpg'));

% Smoother bias plot.
figure;
scatter(fr./dsp, abs(smoother_bias_box), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
scatter(fr./dsp, abs(smoother_bias_gss), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)

legend('box-filtered', 'Gaussian-filtered')
xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
ylabel('$\left|\frac{\delta K}{K}\right|$')
title(strcat('Smoother Bias at $r = $', {' '}, string(fr)))
% saveas(gcf, strcat(img_fdr, 'smoother-bias\bias-', string(fr), 'r', '.jpg'));

% Amplification coefficient plot.
% figure;
% errorbar(flip(resol), flip(amp_box), std(alpha_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
% hold on
% errorbar(flip(resol), flip(amp_gss), std(alpha_gss, 0, 2), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)
% 
% legend('box-filtered', 'Gaussian-filtered')
% xlabel('Normalized Spacing $s$')
% ylabel('$\alpha$')
% title('Quadratic Amplification Coefficient')

% saveas(gcf, strcat(img_fdr, 'amp-coeff\amp-', string(fr), 'r', '.jpg'));
