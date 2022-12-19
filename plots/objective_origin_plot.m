%%%%%%%%%%%%%%%%%% Scatter plot of objective origins %%%%%%%%%%%%%%%%%%%
close all
clear

% Font size for title.

fontName = 'Arial';
% Noise level chosen to be plotted.
noise_idx = 7;

load('objective-momentum.mat');
fsize = 8;

% Unfiltered.
figure;
t = tiledlayout(1,2,'TileSpacing','compact');

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
cbfontsize = 1.5*fsize;

cb = colorbar;
cb.Label.String = '$\frac{|\delta\vec{I}|}{I}$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = cbfontsize;
cb.Label.Rotation = 0;
cb.Label.Units = 'normalized';
cb.Label.Position = [0.75 1.125 0];

xlabel('X/R','fontName',fontName,'FontSize',fsize,'Interpreter','none')
ylabel('Y/R','fontName',fontName,'FontSize',fsize,'Interpreter','none')
zlabel('Z/R','fontName',fontName,'FontSize',fsize,'Interpreter','none')
text(-2.5,2,2.9,'(a)', 'FontSize', fsize,'FontName',fontName,'Interpreter','none','fontWeight','normal')
xlim([-2 2]);
ylim([-2 2]);
zlim([-2 2]);
xticks([-2:1:2]);
yticks([-2:1:2]);
zticks([-2:1:2]);
axis square

nexttile;
dot_size = 20;
% Normalize momentum identity remainder physically.
% scatter3(origin_unf1(:,1), origin_unf1(:,2), origin_unf1(:,3), dot_size, squeeze(res_unf1*vf.fluid.density*vf.scale.len^4/i0), 'filled')
scatter3(origins(:,1), origins(:,2), origins(:,3), dot_size, squeeze(res_unf(noise_idx,:,1)*vf.fluid.density*vf.scale.len^4/i0), 'filled')

cb = colorbar;
cb.Label.String = '$\frac{|\delta\vec{\epsilon}|}{I}$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = cbfontsize;
cb.Label.Rotation = 0;
cb.Label.Units = 'normalized';
cb.Label.Position = [0.75 1.125 0];

xlabel('X/R','fontName',fontName,'FontSize',fsize,'Interpreter','none')
ylabel('Y/R','fontName',fontName,'FontSize',fsize,'Interpreter','none')
zlabel('Z/R','fontName',fontName,'FontSize',fsize,'Interpreter','none')
text(-2.5,2,2.9,'(b)','FontSize', fsize,'FontName',fontName,'Interpreter','none','fontWeight','normal')
xlim([-2 2]);
ylim([-2 2]);
zlim([-2 2]);
xticks([-2:1:2]);
yticks([-2:1:2]);
zticks([-2:1:2]);
axis square


% Omitting resolution error and smoothing biases, which are negligible.
disp('Relative errors with respect to reference origin')
fprintf('Unfilteed: %f\n', err_rel_unf(noise_idx))
fprintf('Box: %f\n', err_rel_box(noise_idx))
fprintf('Gaussian: %f\n', err_rel_gss(noise_idx))

fig = gcf;
fig.Units = 'centimeters';
fig.Position(3) = 17.4;
fig.Position(4) = 7.5;
%exportgraphics(fig,'ObjectiveOrigin.pdf','ContentType','vector','BackgroundColor','None')
%%%%%%%%%%%%%%%% Residual plots %%%%%%%%%%%%%%%%

% Font.
font = 'Arial';
fontSize = 8;

% Plot of mean errors of objective origins.
figure;
t2 = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile
errorbar(props, mean_err_unf, errm_unf_sd, 'ko', 'Color', 'green', 'MarkerFaceColor', 'green', 'LineWidth', 1)
hold on
errorbar(props, mean_err0_unf, err0m_unf_sd, 'd', 'Color', 'blue', 'MarkerFaceColor', 'blue', 'LineWidth', 1)
xlabel('$\frac{\delta u}{u_0}$', 'FontSize', 1.5*fontSize)
ylabel('$\frac{|\delta \vec{I}|}{I}$', 'FontSize', 1.5*fontSize)
legend({'Objective x_o', 'Centroid x_o'},'FontName',font,'FontSize',fontSize,'Interpreter','tex','Location','northwest')
title('(a)', 'FontName', font, 'FontSize', fontSize, 'interpreter', 'none', 'fontweight', 'normal')
ax = gca;
ax.YLabel.Rotation = 0;
ax.YLabel.Position(1) = -0.75;
axis square
xlim([-0.1 3.1])

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
    h1{i}=scatter(mres_unf(n,:), mag_err_unf(n,:), 'filled', cols{i}, 'o');
    hold on
    h2{i}=scatter(mres0_unf(n,:), mag_err0_unf(n,:), 'filled', cols{i}, 'd');
    %h2{i}=scatter(mres0_unf(n,:), mag_err_unf(n,:), 'filled', cols{i}, 'd');
    hold on
end
%hold off
xlabel('$\frac{|\vec{\epsilon}|}{I}$', 'FontSize', 1.5*fontSize)
ylabel('$\frac{|\delta\vec{I}|}{I}$', 'FontSize', 1.5*fontSize)
lgd = repmat(props(noise_idx), 2, 1);

% legend plots and regression lines
m1 = scatter([NaN NaN],[NaN NaN],'filled','k','o');
m2 = scatter([NaN NaN],[NaN NaN],'filled','k','d');

dE = [0:.01:.08];
dIobj = 0.8055*dE + .01483;
dInat = 1.347*dE - 0.005538;
lin1=plot(dE,dIobj,'k');
lin2=plot(dE,dInat,'k--');

%legend(num2cell(strcat('$\delta u=', string(lgd(:)), '$')));
legend([h1{1},h1{2},h1{3},m1,m2,lin1,lin2],{'\deltau = 0.5u_0',...
    '\deltau = 1.5u_0','\deltau= 2.5u_0','Objective x_o','Centroid x_o',...
    'Objective fit','Centroid fit'},...
    'FontName',font,'Interpreter','tex','location','eastoutside');


title('(b)', 'FontName', font, 'FontSize', fontSize, 'interpreter', 'none', 'fontweight', 'normal')
axis equal
ylim([0 0.09])
xlim([0 0.09])

ax = gca;
ax.YLabel.Rotation = 0;
ax.YLabel.Position(1) = -0.02;
box on


fig = gcf;
fig.Units = 'centimeters';
fig.Position(3) = 17.4;
fig.Position(4) = 7;
exportgraphics(fig,'OriginComparison.pdf','ContentType','vector','BackgroundColor','None')