% Generate a panel of 3 figures illustrating Hill's vortex: 3D, xz-contour,
% and 3D with noise.

% Font.
font = 'Arial';
fontSize = 8;
% Only for title.
fontWeight = 'normal';
close all

%%%%%% Plots of 3D Hill's vortex without noise %%%%%%%%
l = 1.5;
vr = 2/3;
% Radius of vortex.
a = l*vr;
spr = 0.25;
u0 = 1;
kappa = a/spr;
[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, u0, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

% Combine figures.
t = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile

% Plot internal and external fields separately.
r = sqrt(sum(vf.X.^2, 4));
Int = repmat(r<=a, 1, 1, 1, 3);
vf.plotVector(vf.U_e.*Int, 0, '', 'Color', 'red');
hold on
vf.plotVector(vf.U_e.*(~Int), 0, '', 'Color', 'cyan')
%title('(a) Hill''s vortex in absolute frame', ...
%    'FontName', font, 'FontSize', fontSize, 'FontWeight', fontWeight, 'interpreter', 'tex')

title('(a)','FontName', font, 'FontSize', fontSize, 'FontWeight', fontWeight, 'interpreter', 'none')

xlabel('x/R', 'FontName', font, 'FontSize', fontSize,'interpreter','none')
ylabel('y/R', 'FontName', font, 'FontSize', fontSize,'interpreter','none')
zlabel('z/R', 'FontName', font, 'FontSize', fontSize,'interpreter','none')
ts = -l:0.5:l;
xticks(ts)
yticks(ts)
zticks(ts)

xlim([-l l])
ylim([-l l])
zlim([-l l])

axis square
axA = gca;
axA.FontName = font;
axA.FontSize = fontSize;
axA.XLabel.Position(3) = axA.XLabel.Position(3)+0.5;
axA.YLabel.Position(3) = axA.YLabel.Position(3)+0.75;

%%%%%%%% xz-cross section %%%%%%%%
% Spacing should divide 1.
%spr = 0.05; messing with spacing for quiver plot
spr = 0.25;
l = 1.5;
vr = 2/3;
sp = spr*l*vr;
[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, 1, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

% Extract xz-plane.
xp = permute(squeeze(x(1,:,:)), [2 1]);
zp = permute(squeeze(z(1,:,:)), [2 1]);
yIdx = floor(size(x,1)/2) + 1;
up = permute(squeeze(u(yIdx,:,:)), [2 1]);
wp = permute(squeeze(w(yIdx,:,:)), [2 1]);

rho = xp.^2 + zp.^2;
Intp = rho<=a;

nexttile
%tp = streamslice(xp,zp,up,wp,2); messing with streamslice settings

tp = streamslice(xp,zp,up,wp,0.5,'noarrows');
hold on
quiver(xp, zp, up.*Intp, wp.*Intp,1, 'Color', 'red','LineWidth',1);
hold on
quiver(xp, zp, up.*(~Intp), wp.*(~Intp),1, 'Color', 'cyan','LineWidth',1);
set(tp, 'LineWidth', 0.25,'Color',[0.5 0.5 0.5])
xticks(ts)
yticks(ts)
xlim([-l l])
ylim([-l l])

xlabel('x/R', 'FontName', font, 'FontSize', fontSize,'interpreter','none')
ylabel('z/R', 'FontName', font, 'FontSize', fontSize,'interpreter','none')
% title('(b) \itxz\rm-cross section of absolute velocity', ...
%     'FontName', font, 'FontSize', fontSize, 'FontWeight', fontWeight, 'interpreter', 'tex')
title('(b)','FontName', font, 'FontSize', fontSize, 'FontWeight', fontWeight, 'interpreter', 'none')
grid off
box on
axis square

axB = gca;
axB.FontName = font;
axB.FontSize = fontSize;

%%%%%% Plot of 3D Hill's vortex with noise %%%%%%%%
n = 1.5;
%n = 1;
N = vf.noise_uniform(n*u0);
Un = vf.U_e + N;
unp = permute(squeeze(Un(yIdx,:,:,1)), [2 1]);
wnp = permute(squeeze(Un(yIdx,:,:,3)), [2 1]);

% Interior to the circular cross section.

nexttile

quiver(xp, zp, unp.*Intp, wnp.*Intp,1, 'Color', 'red','LineWidth',1);
hold on
quiver(xp, zp, unp.*(~Intp), wnp.*(~Intp), 1,'Color', 'cyan','LineWidth',1);

xticks(ts)
yticks(ts)
xlim([-l l])
ylim([-l l])

xlabel('x/R', 'FontName', font, 'fontSize', fontSize,'interpreter','none')
ylabel('z/R', 'FontName', font, 'fontSize', fontSize,'interpreter','none')
%title(sprintf('(c) \\itxz\\rm-cross section with %d%% noise', int32(n*100)), ...
    %'FontName', font, 'FontSize', fontSize, 'FontWeight', fontWeight, 'interpreter', 'tex')
title('(c)','FontName', font, 'FontSize', fontSize, 'FontWeight', fontWeight, 'interpreter', 'none')

axis square

axC = gca;
axC.FontName = font;
axC.FontSize = fontSize;

fig = gcf;
fig.Units = 'centimeters';
fig.Position(3) = 17.4;
fig.Position(4) = 6.5;
exportgraphics(fig,'HillVortex.pdf','ContentType','vector','BackgroundColor','None')
