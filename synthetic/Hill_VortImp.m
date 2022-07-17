function I = Hill_VortImp(vf, R, origin, with_noise)

r = sqrt(sum(vf.X_e.^2, 4));

Int = repmat(r<=R, 1, 1, 1, 3);

if ~isequal(size(origin), [3 1])
    error('Invalid origin')
end

if with_noise
    I = squeeze(vf.fluid.density/2 * ...
        sum(cross(VelocityField.subtract3Vector(vf.X_e, origin), ...
        vf.vorticity(1), 4).*Int, [1 2 3], 'omitnan') * ...
        vf.solver.dv*vf.scale.len);
else
    I = squeeze(vf.fluid.density/2 * ...
        sum(cross(VelocityField.subtract3Vector(vf.X_e, origin), ...
        vf.vort_e, 4).*Int, [1 2 3], 'omitnan') * ...
        vf.solver.dv*vf.scale.len);
end