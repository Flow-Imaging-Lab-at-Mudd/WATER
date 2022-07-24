%%%%%%%%%%%%%%%%%% Scatter plot of objective origins %%%%%%%%%%%%%%%%%%%
% Font size for title.
fsize = 12;
% Noise level chosen to be plotted.
noise_idx = 7;

% Unfiltered.
figure;
t = tiledlayout(1,2);

% There is an OUTLIER in this dataset for 150% noise which is manually
% removed for the presentation of the graph. This will be verbally
% reported.
% origin_ol = squeeze(origin_unf(4, 12, :));
% % New container without this origin.
% origin_unf1 = squeeze(origin_unf(noise_idx,:,:));
% origin_unf1(12,:) = [];
% err_ol = mag_err_unf(4,12);
% mag_err_unf1 = mag_err_unf(4,:);
% mag_err_unf1(12) = [];
% res_ol = res_unf(4,12);
% res_unf1 = squeeze(res_unf(4,:,1));
% res_unf1(12) = [];


nexttile;
dot_size = 20;
origins = squeeze(origin_unf(noise_idx,:,:));
% scatter3(origin_unf1(:,1), origin_unf1(:,2), origin_unf1(:,3), dot_size, mag_err_unf1, 'filled')
scatter3(origins(:,1), origins(:,2), origins(:,3), dot_size, mag_err_unf(noise_idx,:), 'filled')
cbfontsize = 12;

cb = colorbar;
cb.Label.String = '$\frac{|\delta\vec{I}|}{I}$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = cbfontsize;
cb.Label.Rotation = 0;

xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('(a) Error of objective origins', 'FontSize', fsize)
axis square

nexttile;
dot_size = 20;
% Normalize momentum identity remainder physically.
% scatter3(origin_unf1(:,1), origin_unf1(:,2), origin_unf1(:,3), dot_size, squeeze(res_unf1*vf.fluid.density*vf.scale.len^4/i0), 'filled')
scatter3(origins(:,1), origins(:,2), origins(:,3), dot_size, squeeze(res_unf(noise_idx,:,1)*vf.fluid.density*vf.scale.len^4/i0), 'filled')

cb = colorbar;
cb.Label.String = '$\frac{|\delta\vec{\epsilon}|\rho}{I}$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = cbfontsize;
cb.Label.Rotation = 0;

xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('(b) Residual of momentum identity', 'FontSize', fsize)
axis square


% Omitting resolution error and smoothing biases, which are negligible.
disp('Relative errors with respect to reference origin')
fprintf('Unfilteed: %f\n', err_rel_unf(noise_idx))
fprintf('Box: %f\n', err_rel_box(noise_idx))
fprintf('Gaussian: %f\n', err_rel_gss(noise_idx))

%%%%%%%%%%%%%%%% Residual plots %%%%%%%%%%%%%%%%

% Font.
font = 'Arial';
fontSize = 10;

% Plot of mean errors of objective origins.
figure;
t2 = tiledlayout(1,2);

nexttile
errorbar(props, mean_err_unf, errm_unf_sd, 'ko', 'Color', 'green', 'MarkerFaceColor', 'green', 'LineWidth', 1)
hold on
errorbar(props, mean_err0_unf, err0m_unf_sd, 's', 'Color', 'blue', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
xlabel('$\frac{\delta u}{u_0}$', 'FontSize', fontSize)
ylabel('$\frac{|\delta \vec{I}|}{I}$', 'FontSize', fontSize)
legend({'objective origin', 'natural origin'})
title('(a) Impulse errors of objective and natural origins', 'FontName', font, 'FontSize', fontSize, 'interpreter', 'tex', 'fontweight', 'normal')
ax = gca;
ax.YLabel.Rotation = 0;
ax.YLabel.Position(1) = -0.3;
axis square

% Normalize residual.
mres_unf = vf.fluid.density*vf.scale.len^4/i0*squeeze(res_unf(:,:,1));
mres0_unf = vf.fluid.density*vf.scale.len^4/i0*squeeze(res0_unf(:,:,1));

% Compute average residual.
mres_sd = std(mres_unf, 0, 2);
mres_ave = mean(mres_unf, 2);
mres0_sd = std(mres0_unf, 0, 2);
mres0_ave = mean(mres0_unf, 2);

% % Plot of averaged residual over noise levels.
% figure;
% errorbar(props, mres_ave, mres_sd, 'ko', 'MarkerFaceColor', 'red')
% hold on
% errorbar(props, mres0_ave, mres0_sd, 'ko', 'MarkerFaceColor', 'blue')
% xlabel('$\frac{\delta u}{u_0}$')
% ylabel('$\frac{|\delta\vec{\epsilon}|\rho}{I}$')
% legend({'objective origin', 'central origin'})
% title('Residual of momentum identity')
% ax = gca;
% ax.YLabel.Rotation = 0;
% 
% 
% Plot of impulse error and momentum residual.
nexttile
noise_idx = [2, 4, 6];
cols = {'g', 'b', 'r'};

for i = 1: length(noise_idx)
    n = noise_idx(i);
    scatter(mres_unf(n,:), mag_err_unf(n,:), 'filled', cols{i}, 'o')
    hold on
    scatter(mres0_unf(n,:), mag_err_unf(n,:), 'filled', cols{i}, 'd')
    hold on
end
hold off
xlabel('$\frac{|\delta\vec{\epsilon}|\rho}{I}$', 'FontSize', fontSize)
ylabel('$\frac{|\delta\vec{I}|}{I}$', 'FontSize', fontSize)
lgd = repmat(props(noise_idx), 2, 1);
legend(num2cell(strcat('$\delta u=', string(lgd(:)), '$')));
ax = gca;
ax.YLabel.Rotation = 0;
ax.YLabel.Position(1) = -0.01;
title('(b) Momentum identity residual and impulse error', 'FontName', font, 'FontSize', fontSize, 'interpreter', 'tex', 'fontweight', 'normal')
axis square

% fig = gcf;
% fig.Units = 'centimeters';
% fig.Position(3) = 11.9;
% fig.Position(4) = 7;