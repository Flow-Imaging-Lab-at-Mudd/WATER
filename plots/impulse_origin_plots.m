% Load computed data of impulse errors at uniformly sampled impulse
% locations.
close all
load('origin-graph.mat')

% Font.
font = 'Arial';
fontSize = 8;
% Only for title.
fontWeight = 'normal';

t = tiledlayout(1, 2,'TileSpacing','compact','Padding','compact');

% Error magnitude after noise without smoothing.
nexttile;
vfp.plotScalar(mag_err(:,:,:,2), 0, '');
text(-1.5,1,2,'(a)','FontName',font,'FontSize',fontSize,'fontWeight',fontWeight,'Interpreter','none')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = '$\frac{|\delta\vec{I}|}{I}$';
c.Label.Rotation = 0;
c.Label.FontSize = 1.5*fontSize;
xlabel('X/R','FontName',font,'FontSize',fontSize,'Interpreter','none')
ylabel('Y/R','FontName',font,'FontSize',fontSize,'Interpreter','none')
zlabel('Z/R','FontName',font,'FontSize',fontSize,'Interpreter','none')
xticks([-1 -0.5 0 0.5 1]);
yticks([-1 -0.5 0 0.5 1]);

c.Label.Units = 'normalized';
% c.Limits = [0 0.08];

axis square
c.Label.Position = [0.7 1.125 0];

ax = gca;
ax.XTickLabelRotation = 0;

% Plot of standard deviations.
nexttile;
mag_err_sd = squeeze(mag_err_sd);
vfp.plotScalar(mag_err_sd(:,:,:,2), 0, '');
text(-1.5,1,2,'(b)','FontName',font,'FontSize',fontSize,'fontWeight',fontWeight,'Interpreter','none')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = 'SD$\left(\frac{|\delta\vec{I}|}{I}\right)$';
c.Label.Rotation = 0;
c.Label.FontSize = 1.25*fontSize;
xlabel('X/R','FontName',font,'FontSize',fontSize,'Interpreter','none')
ylabel('Y/R','FontName',font,'FontSize',fontSize,'Interpreter','none')
zlabel('Z/R','FontName',font,'FontSize',fontSize,'Interpreter','none')
xticks([-1 -0.5 0 0.5 1]);
yticks([-1 -0.5 0 0.5 1]);
% c.Limits = [0 0.08];
c.Label.Units = 'normalized';
axis square
c.Label.Position = [1.625 1.125 0];
ax1 = gca;
ax1.XTickLabelRotation = 0;

fig = gcf;
fig.Units = 'centimeters';
fig.Position(2) = 5;
fig.Position(3) = 17.4;
fig.Position(4) = 10;
exportgraphics(fig,'OriginSweep.pdf','ContentType','vector','BackgroundColor','None')