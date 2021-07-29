% Verify the expectation that the unfiltered error term in impulse
% computation is independent from the nature of the field.

[x, y, z, u, v, w, Mag] = hill_vortex_3D(0.05, 1, 1, 1);

% Create velocity fields to compare.
vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
vf2 = VelocityField.import_grid_separate(x,y,z,u,v,w);

% Upscale velocity magnitude.
vf.U_e = 1000 * vf.U_e;
vf.deriveQuantities()
vf.updateGlobal()

u1 = vf.meanSpeed(0, 0);
N = vf.noise_uniform(u1);

origin = [0 0.5 0]';

dI1 = vf.impulse(origin, 1) - vf.impulse(origin, 0);

for i = 1: 10
    % Alter velocity of the second field.
    vf2.U_e = 2*u1*rand(size(vf2.U_e)) - u1;
    vf2.deriveQuantities()
    vf2.updateGlobal()
    
    vf2.setNoise(N)
    
    dI2 = vf2.impulse(origin, 1) - vf2.impulse(origin, 0);
    disp(abs((dI1 - dI2) ./ dI1))
end




