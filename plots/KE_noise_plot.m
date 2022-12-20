clear
close all

% Font.
font = 'Arial';
fontSize = 8;

%%%%%% Plots of 3D Hill's vortex without noise %%%%%%%%
l = 1;
vr = 1;
r = l*vr;
% Radius of vortex.
spr = 0.05;
% lower resolution option for supplement.
spr2 = 0.1;

[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, 1, 1);
[x2, y2, z2, u2, v2, w2] = Hill_Vortex(spr2, l, vr, 1, 1);

vf = VelocityField.importCmps(x, y, z, u, v, w);
vf2 = VelocityField.importCmps(x2, y2, z2, u2, v2, w2);

% Proportions of noise.
props = 0: 0.5: 3;

% Compute vortical KE of Hill's vortex.
KEf = @(vf, n) Hill_VortKE(vf, r, n);
K0 = Hill_KE(vf.fluid.density, vf.scale.len, 1, 1);

% Panel figures.
figure;
t = tiledlayout(1,3,'TileSpacing','loose','Padding','compact');

nexttile
% Plot signed error from one run.
KE_err_run(vf, props, K0, KEf, 1, [], {'signed'});
    
box on
axA = gca;
axA.FontName = font;
axA.FontSize = fontSize;
xlim([-0.1 3.1])
ylim([-0.05 0.75])
axA.XLabel.FontSize = 1.5*fontSize;
axA.YLabel.FontSize = 1.5*fontSize;
% axA.YLabel.Rotation = 0;
% axA.YLabel.Position(1) = axA.YLabel.Position(1)-0.025;
title('(a) \kappa = 20','FontName',font,'FontSize',fontSize,'Interpreter','tex')

num_ite = 20;

nexttile
% Plot error magnitude over iterations.
KE_err_run(vf2, props, K0, KEf, num_ite, [], {'mag'});
axB = gca;
axB.FontName = font;
axB.FontSize = fontSize;
xlim([-0.1 3.1])
ylim([0 0.18])
axB.XLabel.FontSize = 1.5*fontSize;
axB.YLabel.FontSize = 1.5*fontSize;
% axB.YLabel.Rotation = 0;
% axB.YLabel.Position(1) = axB.YLabel.Position(1)-0.025;
title('(b) \kappa = 10','FontName',font,'FontSize',fontSize,'Interpreter','tex')

nexttile
% Plot error magnitude over iterations.
KE_err_run(vf, props, K0, KEf, num_ite, [], {'mag'});
axC = gca;
axC.FontName = font;
axC.FontSize = fontSize;
xlim([-0.1 3.1])
ylim([0 0.18])
axC.XLabel.FontSize = 1.5*fontSize;
axC.YLabel.FontSize = 1.5*fontSize;
% axC.YLabel.Rotation = 0;
% axC.YLabel.Position(1) = axC.YLabel.Position(1)-0.025;
title('(c) \kappa = 20','FontName',font,'FontSize',fontSize,'Interpreter','tex')

fig = gcf;
fig.Units = 'centimeters';
fig.Position(3) = 17.4;
fig.Position(4) = 6.5;
exportgraphics(fig,'HillKENoise.pdf','ContentType','vector','BackgroundColor','None')

