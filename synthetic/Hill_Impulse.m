function [I, dIt] = Hill_Impulse(density, len_unit, a, u0, z)
% Analytical result for impulse of a Hill's vortex with a limit of
% integration in y as well as its time derivative.

if ~exist('y', 'var')
    z = a;
end

I = density*pi*15/8*u0/a^2 * ...
    (a^4*z - 2/3*a^2*z.^3 + 1/5*z.^5 + 8/15*a^5)*len_unit^4;

I = [0 0 I]';

dIt = -45/16*density * pi/a^2 * u0^2 * ...
    (a^4* - 2*a^2*z.^2 + z.^4)*len_unit^4;

dIt = [0 0 dIt];

