function [mean_err_box, mean_err_gss, smoother_bias_box, smoother_bias_gss, ...
    amp_box, amp_gss] = KE_feature_radius(sp, radii, props)
% 'sp' is normalized spacing of Hill's Vortex.
% 'fr' is the feature radius of the vortex.

radii_count = length(radii);

% Trials per radius.
num_ite = 5;

% Containers for data across all runs.
err_box = zeros(radii_count, num_ite);
err_gss = zeros(radii_count, num_ite);

bias_box = zeros(radii_count, num_ite);
bias_gss = zeros(radii_count, num_ite);

alpha_box = zeros(radii_count, num_ite);
alpha_gss = zeros(radii_count, num_ite);

for k = 1: radii_count
    % Construct Hill vortex with specified resolution.
    fr = radii(k);
    [x, y, z, u, v, w, ~] = hill_vortex_3D(sp, fr, 1, 1);
    vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
    
    % Focus on the global region with the freestream velocity subtracted.
    vf.addVelocity(-vf.U(1,1,1,:));
    % Zoom into vortical region. Adequate resolution is presumed.
    center = fr * repmat([-1 1], 3, 1);
    % Skip under-resolved case.
    if prod(center(:,2) - center(:,1)) == 0
        break
    else
        vf.setRangePosition(fr * repmat([-1 1], 3, 1));
    end
    
    for i = 1: num_ite    
        % Run script for KE error samrpling.
        [err_box(k, i), err_gss(k, i), bias_box(k, i), ...
            bias_gss(k, i), alpha_box(k, i), alpha_gss(k, i)] = ...
            KE_mean_err_run(vf, props, fr);
    end
end

mean_err_box = mean(err_box, 2);
mean_err_gss = mean(err_gss, 2);

smoother_bias_box = mean(bias_box, 2);
smoother_bias_gss = mean(bias_gss, 2);

amp_box = mean(alpha_box, 2);
amp_gss = mean(alpha_gss, 2);

figure;
errorbar(radii / sp, mean_err_box, std(err_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth',1)
hold on
errorbar(radii / sp, mean_err_gss, std(err_gss, 0, 2), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)

% Mean filter error plot.
legend('box-filtered', 'Gaussian-filtered', 'Interpreter', 'latex')
xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
ylabel('$\bar{\left|\frac{\delta K}{K}\right|}$')
title(strcat('Mean Error of Smoother over $\delta u = $', {' '}, ...
    string(props(1)), '-', string(props(end)*100), '\% $\bar{u}$'))

% Folder for plots.
img_fdr = strcat('C:\Users\derek\flow\trials\ke\feature\s=', string(sp), '\');
mkdir(img_fdr);
saveas(gcf, strcat(img_fdr, 'err-', string(sp), 's', '.jpg'));

% Smoother bias plot.
figure;
scatter(radii / sp, smoother_bias_box, 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
scatter(radii / sp, smoother_bias_gss, 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)

legend('box-filtered', 'Gaussian-filtered', 'Interpreter', 'latex')
xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
ylabel('$\left|\frac{\delta K}{K}\right|$')
title('Smoother Bias')
saveas(gcf, strcat(img_fdr, 'bias-', string(sp), 's', '.jpg'));

% Amplification coefficient plot.
figure;
errorbar(radii / sp, amp_box, std(alpha_box, 0, 2), 'ko', 'MarkerFaceColor','red', 'LineWidth', 1)
hold on
errorbar(radii / sp, amp_gss, std(alpha_gss, 0, 2), 'ko', 'MarkerFaceColor','blue', 'LineWidth', 1)

legend('box-filtered', 'Gaussian-filtered', 'Interpreter', 'latex')
xlabel(strcat('Feature Resolution $\frac{r}{s}$'))
ylabel('$\alpha$')
title('Quadratic Amplification Coefficient')
saveas(gcf, strcat(img_fdr, 'amp-', string(sp), 's', '.jpg'));

close all