% Generate a panel of 3 figures illustrating Hill's vortex: 3D, xz-contour,
% and 3D with noise.

% Font.
font = 'Times New Roman';
fontSize = 10;
% Only for title.
fontWeight = 'normal';

%%%%%% Plots of 3D Hill's vortex without noise %%%%%%%%
l = 1.5;
vr = 2/3;
% Radius of vortex.
a = l*vr;
spr = 0.25;
u0 = 1;

[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, u0, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

% Combine figures.
t = tiledlayout(1,3);
nexttile

% Plot internal and external fields separately.
r = sqrt(sum(vf.X.^2, 4));
Int = repmat(r<=a, 1, 1, 1, 3);
vf.plotVector(vf.U_e.*Int, 0, '', 'Color', 'red');
hold on
vf.plotVector(vf.U_e.*(~Int), 0, '', 'Color', 'cyan')
title('(a) Hill''s vortex in absolute frame', ...
    'FontName', font, 'FontSize', fontSize, 'FontWeight', fontWeight, 'interpreter', 'tex')

xlabel('$x$', 'FontName', font, 'FontSize', fontSize)
ylabel('$y$', 'FontName', font, 'FontSize', fontSize)
ts = -l:0.5:l;
xticks(ts)
yticks(ts)
zticks(ts)

xlim([-l l])
ylim([-l l])
zlim([-l l])

axis square

%%%%%%%% xz-cross section %%%%%%%%
% Spacing should divide 1.
spr = 0.05;
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

nexttile
tp = streamslice(xp,zp,up,wp,2);
set(tp, 'LineWidth', 1.25)
xlim([-l l])
ylim([-l l])

xlabel('$x$', 'FontName', font, 'FontSize', fontSize)
ylabel('$z$', 'FontName', font, 'FontSize', fontSize)
title('(b) \itxz\rm-cross section of absolute velocity', ...
    'FontName', font, 'FontSize', fontSize, 'FontWeight', fontWeight, 'interpreter', 'tex')
grid off

axis square

%%%%%% Plot of 3D Hill's vortex with noise %%%%%%%%
n = 1.5;
N = vf.noise_uniform(n*u0);
Un = vf.U_e + N;
unp = permute(squeeze(Un(yIdx,:,:,1)), [2 1]);
wnp = permute(squeeze(Un(yIdx,:,:,3)), [2 1]);

% Interior to the circular cross section.
rho = xp.^2 + zp.^2;
Intp = rho<=a;
nexttile

quiver(xp, zp, unp.*Intp, wnp.*Intp, 'Color', 'red');
hold on
quiver(xp, zp, unp.*(~Intp), wnp.*(~Intp), 'Color', 'cyan');

xticks(ts)
yticks(ts)
xlim([-l l])
ylim([-l l])

xlabel('$x$', 'FontName', font, 'fontSize', fontSize)
ylabel('$z$', 'FontName', font, 'fontSize', fontSize)
title(sprintf('(c) \\itxz\\rm-cross section with %d%% noise', int32(n*100)), ...
    'FontName', font, 'FontSize', fontSize, 'FontWeight', fontWeight, 'interpreter', 'tex')

axis square
