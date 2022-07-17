function [err] = momentum_origin_obj(origin, vf, thr)
% Objective function for the derivative moment transform of momentum
% 
% Derek Li, July 2021

% Distance units (specified in vf) are not multiplied in this problem.
dv = abs(vf.xsp*vf.ysp*vf.zsp);

% Noisy velocity.
U = vf.U_e + vf.N_e;

% Noisy vorticity.
vort = vf.vorticity(1);

% Threshold vorticity to prevent noise dominance?
vormag = sqrt(sum(vort.^2,4));
vort(:,:,:,1,:) = vort(:,:,:,1,:).*(vormag > thr);
vort(:,:,:,2,:) = vort(:,:,:,2,:).*(vormag > thr);
vort(:,:,:,3,:) = vort(:,:,:,3,:).*(vormag > thr);

% Relative position to the given origin.
X_rel = VelocityField.operate3Vector(vf.X_e, origin, @minus);

% Remainder of the momentum identity, norm to be
% minimized.
momentum = dv*squeeze(sum(U, [1 2 3], 'omitnan'));
impulse = 0.5*dv*squeeze(sum(cross(X_rel, vort, 4), [1 2 3], 'omitnan'));
flux = 0.5*vf.intCubicSurf_doublecross(X_rel, U);

rem1 = momentum - impulse + flux;

r1 = norm(rem1);

err = r1;

% % Compute gradient for this objective function.
% derv = zeros(3, 1);
% % Identity matrix used for diff indexing.
% I = eye(3);

% for i = 1: 3
%     derv1 = rem1'/r1 * cross(I(:,i), squeeze(dv*sum(vort, [1 2 3], 'omitnan')) + ...
%         vf.intCubicSurf_cross(U));
%     derv(i) = derv1;
% end
