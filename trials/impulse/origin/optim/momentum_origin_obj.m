function [err, derv] = momentum_origin_obj(origin, vf)
% Objective function for the two integral equations that determine the
% Ringuette objective origin, given in (De Voria, 2014).
% 
% Derek Li, July 2021

% Distance units (specified in vf) are not multiplied in this problem.
dv = abs(vf.xsp*vf.ysp*vf.zsp);

% Noisy velocity.
U = vf.U_e + vf.N_e;
% Noisy velocity gradient.
ugrad = vf.gradient(U);
% Noisy vorticity.
vort = vf.vorticity(1);
% Relative position to the given origin.
X_rel = VelocityField.operate3Vector(vf.X_e, origin, @minus);

% Remainder of the momentum identity, norm to be
% minimized.
rem1 = 2*dv*squeeze(sum(U, [1 2 3], 'omitnan')) - ...
    dv*squeeze(sum(cross(X_rel, vort, 4), [1 2 3], 'omitnan')) + ...
    vf.intCubicSurf_doublecross(X_rel, U);

r1 = norm(rem1);

err = r1;

% Compute gradient for this objective function.
derv = zeros(3, 1);
% Identity matrix used for diff indexing.
I = eye(3);

for i = 1: 3
    derv1 = rem1'/r1 * cross(I(:,i), squeeze(dv*sum(vort, [1 2 3], 'omitnan')) + ...
        vf.intCubicSurf_cross(U));
    derv(i) = derv1;
end
