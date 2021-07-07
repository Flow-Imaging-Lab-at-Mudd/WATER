% Paremeters held constant.
sp = 0.1;
fr = 1;
u0 = 1;

% Proportional noise.
props = 0: 0.1: 3;

[x, y, z, u, v, w, ~] = hill_vortex_3D(sp, fr, u0, 1);

vf = VelocityField.import_grid_separate(x,y,z,u,v,w);

% Sample origins uniformly from the grid. This is done by making another
% velocity field, which has no velocity content. The positions of this
% velocity field are the origins.
osp = 1 * ones(1, 3);
oends = [-2 2; -2 2; -2 2];

clear X
[X(:,:,:,1), X(:,:,:,2), X(:,:,:,3)] = meshgrid(oends(1,1): osp(1): oends(1,2), ...
    oends(2,1): osp(2): oends(2,2), oends(3,1): osp(3): oends(3,2));
vfp = VelocityField(X, zeros(size(X)));

% Container for optimal solutions found upon different initial origin
% guesses.

% Initial guesses of origin.
[o, min_err] = fminsearch(@(x) min_err_origin(x, 'box', vf, props, fr, u0), [2 2 2]')