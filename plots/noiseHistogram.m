% histogram of error due to noise

clear
close all

% Constant parameters.
l = 1;
vr = 1;
fr = l*vr;
u0 = 1;
min_fres = 10;
spr = 1 / min_fres;
u_mean = u0;

% Windowing paramters for 2-parameter noise model.
beta = 0.25;
win = 16;
op = 0.75;

% Noise levels.
props = [1.5];

% Vortex parameters.
density = 1000;
len_unit = 1e-3;
% Theoretical impulse values.
I0 = Hill_Impulse(density, len_unit, fr, u0);
I0mag = norm(I0);

[x, y, z, u, v, w] = Hill_Vortex(spr, l, vr, u0, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);
% N = vf.noise_uniform(props(1)*u_mean);
N = vf.noise_localcor(1, win, op, beta);

origin = [0 0 0]';

I0vox = vf.fluid.density/2 * ...
                    cross(VelocityField.subtract3Vector(vf.X_e, origin), ...
                        vf.vorticity(1), 4) * ...
                    vf.solver.dv*vf.scale.len;

INvox = vf.fluid.density/2 * ...
                    cross(VelocityField.subtract3Vector(vf.X_e, origin), ...
                        vf.vort_e, 4) * ...
                    vf.solver.dv*vf.scale.len;

err = (INvox - I0vox)/I0mag;

bins = [-0.0013:.0001:0.0013];

figure
t = tiledlayout(1,3,'tilespacing','compact','padding','compact')

nexttile
histogram(err(:,:,:,1),bins)
title('X','fontName','Arial','fontSize',8,'FontWeight','normal','Interpreter','none')
xlabel('Voxel Error (normalized)','fontName','Arial','fontSize',8,'FontWeight','normal','Interpreter','none')

axX = gca;
axX.XAxis.Exponent = 0;

nexttile
histogram(err(:,:,:,2),bins)
title('Y','fontName','Arial','fontSize',8,'FontWeight','normal','Interpreter','none')
xlabel('Voxel Error (normalized)','fontName','Arial','fontSize',8,'FontWeight','normal','Interpreter','none')

axY = gca;
axY.XAxis.Exponent = 0;

nexttile
histogram(err(:,:,:,3),bins)
title('Z','fontName','Arial','fontSize',8,'FontWeight','normal','Interpreter','none')
xlabel('Voxel Error (normalized)','fontName','Arial','fontSize',8,'FontWeight','normal','Interpreter','none')

axZ = gca;
axZ.XAxis.Exponent = 0;

fig = gcf
fig.Units = 'centimeters';
fig.Position(3) = 14;
fig.Position(4) = 6;
exportgraphics(fig,'ImpulseNoiseHistogram.pdf','ContentType','vector','BackgroundColor','None')
