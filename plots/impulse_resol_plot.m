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
err_level = 0.1;
% Iterations.
num_ite = 10;

% Resolutions.
min_fres = 20;
fres_inc = 20;
max_fres = 160;

% Vortex parameters.
density = 1000;
len_unit = 1e-3;

% Theoretical values.
origin = [0 0 0]';
I0 = Hill_Impulse(density, len_unit, r, u0, r);

[fres, dI, dI_box, dI_gss, dI0, Ibias_box, Ibias_gss, di, di_box, ...
    di_gss, di0, mag_bias_box, mag_bias_gss, dI_sd, dI_box_sd, dI_gss_sd,...
    di_sd, di_box_sd, di_gss_sd, vfds, t] = ...
    impulse_resol(l, vr, u0, min_fres, max_fres, fres_inc, origin, props, err_level, num_ite, window_params, true);

% % Organize into tiles.
% figure;
% t = tiledlayout(1,2);
% axes{1}.Parent = t;
% axes{1}.Layout.Tile = 1;
% % axes{1}.Legend = legend({'unfiltered', 'box-filtered', 'Gaussian-filtered'})
% 
% axes{2}.Parent = t;
% axes{2}.Layout.Tile = 2;
% % legend(axes{2}, {'unfiltered', 'box-filtered'})

% [~, dK, dK_box, dK_gss, Kbias_box, Kbias_gss, dK0, dK_sd, dK_box_sd, dK_gss_sd] = ...
%     KE_resol(l, vr, u0, min_fres, max_fres, fres_inc, props, err_level, num_ite, window_params, true, vfds);