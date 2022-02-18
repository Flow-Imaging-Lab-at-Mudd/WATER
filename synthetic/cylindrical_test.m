% Radial spacing.
rsp = 0.05;
% Angular spacing.
tsp = 0.1;
% z spacing.
zsp = 0.05;

% Construct cylindrical representation.
[r, theta, z, ur, ut, uz] = Hill_Vortex_Cyl(rsp, tsp, zsp, 1, 1, 0);
vf = CylVF.importCmps(r, theta, z, ur, ut, uz, 0);

% Compute KE.
disp('Theoretical KE:')
K = Hill_KE(vf.fluid.density, vf.scale.len, 1, 1);
disp(K)
disp('Computed KE:')
disp(vf.kineticEnergy())

% Plot velocity field by converting to Cartesian.
x = r.*cos(theta);
y = r.*sin(theta);
u = ur.*cos(theta);
v = ur.*sin(theta);
w = uz;
vf = VelocityField.importCmps(x, y, z, u, v, w, 0);

vf.plotVector(vf.U_e, 0, '$\vec{u}$');