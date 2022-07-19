% Comparison of error profiles of impulse and KE computations with respect
% to windowing size.

clear all
close all
startup;

% Constant parameters.
%l = 1.125;
% max field size to keep to a 512^3 voxel fft is 448, leave one window
% outside vortex for all trials
l = 447/320;
%vr = 8/9;
%vr = 6/9;
vr = 320/447;
r = l*vr;
u0 = 1;
spr = 1/160;

kappa = r/spr;

% Number of vectors in the end which windowing does not sample from.
drem = @(n,w,o) mod((n-w),round((1-o)*w));
 
[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, u0, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w, 1);

% Theoretical values of impulse & KE.
I0 = Hill_Impulse(vf.fluid.density, vf.scale.len, r, u0);
% % Compute only vortical impulse within sphere.
% Ih = @(vf, with_noise) Hill_VortImp(vf, 1, [0 0 0]', with_noise);
% Origin used for impulse calculation.
origin = [0 0 0]';
% Windowing sizes.
winsizes = [16 32 48 64];

% Windowing parameters; window sizes for 0.5 overlap.
% winsizes = [16 32 54 64];
overlap = 0.5;
% Show the dimensions of downsampled field, assuming the field has uniform
% dimensions.
downsampled_dim(vf.dims(1), winsizes, overlap);
kappaWin = r./(winsizes*spr*(1-overlap));
kappaWin25 = r./(winsizes*spr*(1-0.25));
kappaWin75 = r./(winsizes*spr*(1-0.75));

% Noise level.
props = [0 1.5];
%props = 0;

% Panel figures.
figure;
t = tiledlayout(1,2);

font = 'Arial';
fontSize = 8;

% Compute baseline resolution errors.
nexttile
impulse_winsize(vf, I0, origin, props, winsizes, overlap, 1, {'dim', 'resol'});
title('(a)','fontName',font,'fontSize',fontSize,'interpreter','none','fontWeight','normal')
legend({'Unfiltered', 'Gaussian'},'fontName',font,'fontSize',fontSize-1,'interpreter','none','location','southwest')
% Empirically fixed.

box on
% % Compute error magnitude under noise.
% nexttile
% impulse_winsize(vf, I0, origin, props, winsizes, overlap, 20, {'mag'});
% title('(b) Impulse error magnitude under noise')
% ylim([0 0.17])
% box on
xlim([winsizes(1)-2 winsizes(end)+2])
ylim([-0.16 0.02])
% Compute windowing errors of different overlaps.
nexttile;
% o = 0.25.
% winsizes1 = [16 32 50 73];
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ax] = impulse_winsize(vf, I0, origin, props, winsizes, 0.25, 1, {'dim', 'win'});
ax{1}.Marker = 's';
ax{1}.MarkerFaceColor = 'none';
ax{1}.MarkerEdgeColor = 'black';
hold on
% o = 0.5.
% winsizes2 = [16 32 54 64];
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ax] = impulse_winsize(vf, I0, origin, props, winsizes, 0.5, 1, {'dim', 'win'});
ax{1}.Marker = 'd';
ax{1}.MarkerFaceColor = 'none';
ax{1}.MarkerEdgeColor = 'black';
hold on
% o = 0.75.
% winsizes3 = [16 32 42 64];
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ax] = impulse_winsize(vf, I0, origin, props, winsizes, 0.75, 1, {'dim', 'win'});
ax{1}.Marker = 'v';
ax{1}.MarkerFaceColor = 'none';
ax{1}.MarkerEdgeColor = 'black';
title('(b)','fontName',font,'fontSize',fontSize,'interpreter','none','fontWeight','normal')
legend({'o=0.25', 'o=0.5', 'o=0.75'},'fontName',font,'fontSize',fontSize-1,'interpreter','none','location','eastoutside')
ylim([-0.16 0.02])
box on
xticks(winsizes)
xlim([winsizes(1)-2 winsizes(end)+2])
% xticks(unique([winsizes1 winsizes2 winsizes3]))

fig = gcf;
fig.Units = 'centimeters';
fig.Position(3) = 11.9;
fig.Position(4) = 5.5;
exportgraphics(fig,'HillImpulseWin.pdf','ContentType','vector','BackgroundColor','None')
