function g = Hill_VGL2(u0, a, X)
% Given a 4D array in the usual format, compute the L2 norm of the vector
% gradient of Hill's vortex on such positions.

% By the vortex's cylindrical symmetry, the only dependence is on z^2 and r^2.
r2 = sum(X.^2, 4);
z2 = X(:,:,:,3).^2;

% Gradient norm sqaured within the sphere.
g = 9/4*u0^2/a^4*(17*r2 + 6*z2);
% Reset gradient norm squared exterior to and on the sphere.
Ext = (r2 > a^2);
Sph = (r2 == a^2);

Ge = 9/2*u0^2*a^6*(r2 + 3*z2)./(r2 + z2).^5;
g(Ext) = Ge(Ext);

% The gradient is discontinuous at the spherical boundary, so we use the
% average.
Gs = 9/8*u0^2/a^4*(12*a^2 + 7*r2);
g(Sph) = Gs(Sph);

g = sqrt(g);






