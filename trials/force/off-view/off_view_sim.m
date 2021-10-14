% Simulation of a Hill's vortex moving out of view. Impulse and derivatives
% are computed and compared to the theoretical.

sp = 0.05;
r0 = 1;
u0 = 1;

[x, y, z, u, v, w, Mag] = Hill_Vortex(sp, r0, u0, 1);
vf = VelocityField.importCmps(x,y,z,u,v,w);
origin = [0 0 0]';

dt = 0.05;
uc = 3/2 * u0;
T = 0: dt: 2*r0/uc+dt;

% Initially in view.
y = 1;
I_meas = zeros(3, size(T, 2));
I_meas(:,1) = vf.impulse(origin, 0);

for i = 2: size(T, 2)
    t = T(i);
    y = motion_off_view(vf, u0, dt, y);
    I_meas(:,i) = vf.impulse(origin, 0);
end

Iy = I_meas(2,:);

% Compute theoretical values.
Y = 1 - uc*T;
[I_theo, dIt_theo] = theoImpulse(vf, r0, u0, Y);

figure;
plot(T, I_theo)
hold on
scatter(T, Iy, 'r', 'filled')
xlabel('Time ')
ylabel('$I_y$')
title(sprintf('Moving off view at $r/s=%f$', 1/fr))
legend({'theoretical', 'computed'})

% Excerpting the near 0 values of impule in the end.
term = int32(size(Iy, 2)*0.7);
Iye = Iy(1: term);
I_theoe = I_theo(1: term);

figure;
scatter(Y(1:term), abs(Iye-I_theoe)./I_theoe, 'r', 'filled')
xlabel('Axial Terminus')
ylabel('$\delta I_y$')
title('Error percentage of computed impulse')

% Computed derivative.
dIy_dt = vf.scalarTimeDeriv(Iy, dt, 1, 2, 0);
figure;
plot(T, dIt_theo)
hold on
scatter(T, dIy_dt, 'b', 'filled')
xlabel('Time')
ylabel('$\dot{I}_y$')
title('Time derivative of impulse')