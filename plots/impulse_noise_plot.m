% Font.
font = 'Times New Roman';
fontSize = 10;
% Only for title.
fontWeight = 'normal';

%%%%%% Plots of 3D Hill's vortex without noise %%%%%%%%
l = 1;
vr = 1;
% Radius of vortex.
spr = 0.05;

[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, 1, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

% Proportions of noise.
props = 0: 0.5: 6;

I0 = Hill_Impulse(vf.fluid.density, vf.scale.len, 1, 1);

% Panel figures.
t = tiledlayout(1,2);

nexttile
% Plot signed error from one run.
impulse_err_run(vf, props, [0 0 0]', I0, 1, [], {'dim'});
title('(a) Impulse error in $\hat{z}$')
box on

nexttile
num_ite = 20;
% Plot error magnitude over iterations.
impulse_err_run(vf, props, [0 0 0]', I0, 20, [], {'mag'});
title('(b) Magnitude of impulse error')
