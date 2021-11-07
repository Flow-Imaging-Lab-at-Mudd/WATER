function [I, dIt] = theoImpulse(vf, r0, u0, y)
% Analytical expression for impulse of a Hill's vortex with a limit of
% integration in y as well as its time derivative.

if ~exist('y', 'var')
    y = 1;
end

I = vf.fluid.density*pi*15/8*u0/r0^2 * ...
    (r0^4*y - 2/3*r0^2*y.^3 + 1/5*y.^5 + 8/15*r0^5)*vf.scale.len^4;

I = [0 I 0]';

dIt = -45/16* vf.fluid.density * pi/r0^2 * u0^2 * ...
    (r0^4* - 2*r0^2*y.^2 + y.^4)*vf.scale.len^4;

dIt = [0 dIt 0];

