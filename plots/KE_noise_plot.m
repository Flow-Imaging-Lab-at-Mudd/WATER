clear

% Font.
font = 'Times New Roman';
fontSize = 10;

%%%%%% Plots of 3D Hill's vortex without noise %%%%%%%%
l = 1;
vr = 1;
r = l*vr;
% Radius of vortex.
spr = 0.05;

[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, 1, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

% Proportions of noise.
props = 0: 0.5: 3;

% Compute vortical KE of Hill's vortex.
KEf = @(vf, n) Hill_VortKE(vf, r, n);
K0 = Hill_KE(vf.fluid.density, vf.scale.len, 1, 1);

% Panel figures.
figure;
t = tiledlayout(1,2);

nexttile
% Plot signed error from one run.
KE_err_run(vf, props, K0, KEf, 1, [], {'signed'});
    
% title('Signed kinetic energy error', 'FontName', font, 'FontSize', fontSize)
box on

nexttile
num_ite = 10;
% Plot error magnitude over iterations.
KE_err_run(vf, props, K0, KEf, 20, [], {'mag'});
