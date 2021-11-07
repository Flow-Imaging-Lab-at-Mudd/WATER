function [X, Y, Z, u, v, w, Mag] = Hill_Vortex(sp, a, us, z_factor, rem)

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
% us: freestream velocity outside the vortex ring region (in y direction,
% same direction as vortex ring axis), must be nonzero
% z_factor: scale factor on z-direction vector spacing; set to 1 for uniform spacing in x, y, and z
% rem: remove free stream velocity.

% Outputs:
% X,Y,Z: Arrays of X,Y,Z coordinates in meshgrid format. Vortex ring axis
% is in y direction
% u,v,w: Arrays of velocity components in meshgrid format
% Mag: Array of velocity magnitudes

% range of values for x,y,z
% permits different scaling of z-direction resolution using z_factor input
x = [-1:sp:1];
y = [-1:sp:1];
z = [-1:sp*z_factor:1];
[X,Y,Z]=meshgrid(x,y,z);

% default vortex ring center location at center of volume
xc = 0;
yc = 0;
zc = 0;

% get radial distance R to center point for all coordinates
R = sqrt((X-xc).^2 + (Y-yc).^2 + (Z-zc).^2);

% distance from center axis (simplifies coordinates because flow is
% axisymmetric)
r = sqrt((X-xc).^2 + (Z-zc).^2);
theta = atan2(Z,X);

% Inner velocity field (axisymmetric)
A = 15/2*us/a^2;

V = -A/10*(4*r.^2 + 2*Y.^2 - 2*a^2);
U = A*r.*Y/5;

% velocity magnitude
Mag = sqrt(V.^2 + U.^2);

% outer region flow field matches flow over a sphere
Vo = us*((a^2./(r.^2 + Y.^2)).^(5/2).*(2*Y.^2 - r.^2)/(2*a^2)-1);
Uo = 3/2*us/a^2*r.*Y.*(a^2./(r.^2 + Y.^2)).^(5/2);

% outer region flow field is only valid in outer region
Vo (R <= a) = NaN;
Uo (R <= a) = NaN;

% velocity magnitude in outer region
Mo = sqrt(Vo.^2 + Uo.^2);

% Use outer region flow field outside vortex ring radius
V(R>a)=Vo(R>a);
U(R>a)=Uo(R>a);

% Coordinate conversion
u = U.*cos(theta);
v = V;
w = U.*sin(theta);

if exist('rem', 'var') && rem == 1
    v = v + us;
end
