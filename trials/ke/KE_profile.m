function K = KE_profile(vf, with_noise)
% Element-wise kinetic energy on grid.

    K = 1/2*vf.fluid.density * vf.solver.dv * vf.scale.len^2 * ...
                            sum((vf.U_e + with_noise*vf.N_e).^2, 4);