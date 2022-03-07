function Kv = Hill_VortKE(vf, R, with_noise)

r = sqrt(sum(vf.X.^2, 4));

Int = repmat(r<=R, 1, 1, 1, 3);

Kv = 1/2*vf.fluid.density * abs(vf.solver.dv) * vf.scale.len^2 * ...
    sum(((vf.U_e + with_noise*vf.N_e).*Int).^2, 'all', 'omitnan');