function [r, theta, z, ur, ut, uz] = Hill_Vortex_Cyl(rsp, thetasp, zsp, a, u0, rem)

if a > 1
    error('Expected a vortical radius < 1!')
end

% Parameterized in cylindrical coordinates.
r = 0: rsp: 1;
theta = 0: thetasp: 2*pi;
z = -1: zsp: 1;
[r, theta, z] = meshgrid(r, theta, z);

% get radial distance R to center point for all coordinates
R = sqrt(r.^2 + z.^2);

% Interior velocityt field.
% Velocity in z direction.
uz = 3/2*u0*(1 - (2*r.^2 + z.^2)/a^2);
% Velocity in radial direction.
ur = 3/2*u0/a^2*r.*z;
% Tangential velocity is zero.
ut = zeros(size(ur));

% Exterior velocity field.
Uz = u0*((a^2./(z.^2+r.^2)).^(5/2) .* (2*z.^2-r.^2)/(2*a^2) - 1);
Ur = 3/2*u0/a^2*z.*r .* (a^2./(z.^2 + r.^2)).^(5/2);

% Combine interior and exterior, assuming convergence at the spherical
% boundary.
ext = R>a;
uz(ext) = Uz(ext);
ur(ext) = Ur(ext);

if exist('rem', 'var') && rem == 1
    uz = uz + u0;
end
