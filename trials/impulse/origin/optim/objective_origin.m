function err = objective_origin(origin, vf)
% Objective function for the two integral equations that determine the
% Ringuette objective origin.

% Distance units (specified in vf) are not multiplied in this problem.
dv = abs(vf.xresol*vf.yresol*vf.zresol);

% Noisy velocity.
U = vf.U_e + vf.N_e;
% Noisy velocity gradient.
ugrad = vf.gradient(U);
% Noisy vorticity.
vort = vf.vorticity(1);
% Relative position to the given origin.
X_rel = operate3Vector(vf.X_e, origin, @minus);
% Compute the matrix product of velocity gradient and velocity.
dux = NaN([vf.span 3]);

for j = 1: size(vf.X_e, 1)
    for i = 1: size(vf.X_e, 2)
        for k = 1: size(vf.X_e, 3)
            dux(j,i,k,:) = squeeze(ugrad(j,i,k,:,:)) * squeeze(U(j,i,k,:));
        end
    end
end

% Remainders of the two equations, both vectors, their norms to be
% minimized.
rem1 = 2*dv*squeeze(sum(U, [1 2 3])) - ...
    dv*squeeze(sum(cross(X_rel, vort, 4), [1 2 3])) - ...
    vf.intCubicSurf_cross(cross(X_rel, U, 4));

rem2 = -vf.intCubicSurf_vec(U.^2) + ...
    vf.intCubicSurf_cross(cross(X_rel, dux + cross(U, vort, 4), 4));

err = (norm(rem1) + norm(rem2)) / 2;
