function vf = Hill_Plane(sp, a, us, rem)
% A cross section of the Hill Vortex at the xy plane. This is merely a
% variation of the original Hill Vortex construction file, with
% simplification of z_factor = 1 and returning a field object instead of
% position and velocity components.

% range of values for x,y,z
% permits different scaling of z-direction resolution using z_factor input
x = -1: sp: 1;
y = -1: sp: 1;
z = 0;
[X, Y, Z] = meshgrid(x,y,z);

% default vortex ring center location at center of volume
xc = 0;
yc = 0;
zc = 0;

% get radial distance R to center point for all coordinates
R = sqrt((X-xc).^2 + (Y-yc).^2);

% distance from center axis (simplifies coordinates because flow is
% axisymmetric)
r = sqrt((X-xc).^2);
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

vf = VelocityField.importCmps(X, Y, Z, u, v, w);

