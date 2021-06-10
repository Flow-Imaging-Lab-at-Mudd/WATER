vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
vf.setRange([27 30; 20 22; 10 12]);
vf.data.speed = sqrt(sum(vf.U.^2, 4));
% Max speed in region of interest.
speed_r = vf.subsetVector(vf.data.speed);
max_speed_r = max(speed_r, [], 'all');

% Introduce noise proportionally.
props = (1:30) * 0.1;

dK = zeros(size(props));
% Kinetic energy without noise.
k = vf.kineticEnergy(0);

for i = 1: size(props, 2)
    vf.clearNoise();
    vf.noise_uniform(props(i)*361);
    dK(i) = vf.kineticEnergy(1) - k;
end

scatter(props, dK)
xlabel('$\delta u$ (per maximum speed)')
ylabel('$\delta K$ (J)')