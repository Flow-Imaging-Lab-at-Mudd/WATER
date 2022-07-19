% Comparison of error profiles of impulse and KE computations with respect
% to windowing size.

clear
close all
startup;

% Constant parameters.
%l = 1.125;
l = 447/320;
vr = 320/447;
r = l*vr;
u0 = 1;
spr = 1/160;

kappa = r/spr;

% Number of vectors in the end which windowing does not sample from.
drem = @(n,w,o) mod((n-w),round((1-o)*w));
 
[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, u0, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

% Theoretical values of impulse & KE.
K0 = Hill_KE(vf.fluid.density, vf.scale.len, r, u0);
% Compute only vortical impulse.
KEf = @(vf, with_noise) Hill_VortKE(vf, r, with_noise);
% Windowing sizes.
winsizes = [16 32 48 64];

% Windowing parameters; window sizes for 0.5 overlap.
% winsizes = [16 32 54 64];
overlap = 0.5;
% Show the dimensions of downsampled field, assuming the field has uniform
% dimensions.
downsampled_dim(vf.dims(1), winsizes, overlap);

% Noise level.
props = [0 1.5];

% Panel figures.
figure;
t = tiledlayout(1,2);

font = 'Arial';
fontSize = 8;

% Compute baseline resolution errors.
nexttile
KE_winsize(vf, K0, KEf, props, winsizes, overlap, 1, {'resol'});
title('(a)','fontName',font,'fontSize',fontSize,'interpreter','none','fontWeight','normal')
% Empirically fixed.
ylim([-0.15 0])
box on
xlim([winsizes(1)-2 winsizes(end)+2])
% % Compute error magnitude under noise.
% nexttile
% impulse_winsize(vf, I0, origin, props, winsizes, overlap, 10, {'mag'}, KEf);
% title('(b) Impulse error magnitude under noise')
%ylim([-0.15 0])
% box on

% Compute windowing errors of different overlaps.
nexttile;
% o = 0.25.
% winsizes1 = [16 32 50 73];
[~, ~, ~, ~, ~, ~, ~, ~, ~, ax] = KE_winsize(vf, K0, KEf, props, winsizes, 0.25, 1, {'win'});
ax{1}.Marker = 's';
ax{1}.MarkerFaceColor = 'none';
ax{1}.MarkerEdgeColor = 'black';
hold on
% o = 0.5.
% winsizes2 = [16 32 54 64];
[~, ~, ~, ~, ~, ~, ~, ~, ~, ax] = KE_winsize(vf, K0, KEf, props, winsizes, 0.5, 1, {'win'});
ax{1}.Marker = 'd';
ax{1}.MarkerFaceColor = 'none';
ax{1}.MarkerEdgeColor = 'black';
hold on
% o = 0.75.
% winsizes3 = [16 32 42 64];
[~, ~, ~, ~, ~, ~, ~, ~, ~, ax] = KE_winsize(vf, K0, KEf, props, winsizes, 0.75, 1, {'win'});
ax{1}.Marker = 'v';
ax{1}.MarkerFaceColor = 'none';
ax{1}.MarkerEdgeColor = 'black';
title('(b)','fontName',font,'fontSize',fontSize,'interpreter','none','fontWeight','normal')
legend({'o=0.25', 'o=0.5', 'o=0.75'},'fontName',font,'fontSize',fontSize-1,'interpreter','none','location','southwest')

% ylim([-0.17 0])
box on
xticks(winsizes)
xlim([winsizes(1)-2 winsizes(end)+2])

fig = gcf;
fig.Units = 'centimeters';
fig.Position(3) = 11.9;
fig.Position(4) = 6.5;
exportgraphics(fig,'HillKEWin.pdf','ContentType','vector','BackgroundColor','None')

% xticks(unique([winsizes1 winsizes2 winsizes3]))
