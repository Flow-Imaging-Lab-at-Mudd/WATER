function [x, y, z, u, v, w] = Hill_Vortex(sp, a, u0, rem)
% This function outputs synthetic 3D vortex ring flow fields at a
% user-specified resolution
% Vortex rings generated have a velocity field matching that of Hill's
% spherical vortex
% (https://royalsocietypublishing.org/doi/pdf/10.1098/rsta.1894.0006)

% Written by Leah Mendelson 5/1/14
% Last updated: 6/11/21

% Inputs:
% sp: spacing between vectors (normalized by vortex ring radius)
% a: outer radius of sphere / outer radius of vortex region
% u0: freestream velocity outside the vortex ring region (in y direction,
% same direction as vortex ring axis), must be nonzero
% rem: remove free stream velocity.

% Outputs:
% X,Y,Z: Arrays of X,Y,Z coordinates in meshgrid format. Vortex ring axis
% is in y direction
% u,v,w: Arrays of velocity components in meshgrid format
% Mag: Array of velocity magnitudes

if a > 1
    error('Expected a vortical radius < 1!')
end

% range of values for x,y,z
% permits different scaling of z-direction resolution using z_factor input
x = -1: sp: 1;
y = -1: sp: 1;
z = -1: sp: 1;
[x, y, z] = meshgrid(x, y, z);

% default vortex ring center location at center of volume
xc = 0;
yc = 0;
zc = 0;

% get radial distance R to center point for all coordinates
R = sqrt((x-xc).^2 + (y-yc).^2 + (z-zc).^2);

% distance from central z axis (simplifies coordinates because flow is
% axisymmetric)
r = sqrt((x-xc).^2 + (y-zc).^2);
theta = atan2(y, x);

% Interior velocityt field.
% Velocity in z direction.
w = 3/2*u0*(1 - (2*r.^2 + z.^2)/a^2);
% Velocity in planar radial direction.
p = 3/2*u0/a^2*r.*z;

% Exterior velocity field.
wo = u0*((a^2./(z.^2+r.^2)).^(5/2) .* (2*z.^2-r.^2)/(2*a^2) - 1);
po = 3/2*u0/a^2*z.*r .* (a^2./(z.^2 + r.^2)).^(5/2);

% Combine interior and exterior, assuming convergence at the spherical
% boundary.
ext = R>a;
w(ext) = wo(ext);
p(ext) = po(ext);

u = p.*cos(theta);
v = p.*sin(theta);

if exist('rem', 'var') && rem == 1
    v = v - u0;
end
