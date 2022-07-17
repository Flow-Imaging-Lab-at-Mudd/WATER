clear

% Constant parameters.
l = 1;
vr = 1;
r = l*vr;
u0 = 1;

% Windowing parameters.
% winsize = 16;
% overlap = 0.5;
% window_params = [winsize overlap];
window_params = [];

% Noise levels.
props = [0 1.5];
% Desired level of error.
err_level = 0.05;
% Iterations.
num_ite = 10;

% Resolutions.
min_fres = 4;
fres_inc = 4;
max_fres = 40;

% Vortex parameters.
density = 1000;
len_unit = 1e-3;

% Theoretical values.
origin = [0 0 0]';
I0 = Hill_Impulse(density, len_unit, r, u0, r);

% Make tiled figure.
figure;
t = tiledlayout(1,3);

% Resolution error plot.
nexttile;
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, vfds] = KE_resol(l, vr, u0, min_fres, max_fres, fres_inc, props, err_level, num_ite, window_params, {'bias'});

% Magnitude plot.
nexttile;
KE_resol(l, vr, u0, min_fres, max_fres, fres_inc, props, err_level, num_ite, window_params, {'mag'}, vfds);

% Mean signed error plot for high resolutions.
nexttile;
min_fres = 20;
fres_inc = 20;
max_fres = 160;
KE_resol(l, vr, u0, min_fres, max_fres, fres_inc, props, err_level, num_ite, window_params, {'signed'});
