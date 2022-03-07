function [Kv, Ke, K] = Hill_KE(density, len_unit, a, u0)
% Computes the full kinetic energy of a Hill's vortex with the given
% specifications, in the order of vortical, external, and total.

Kv = 23/21*pi*density.*a.^3.*u0.^2.*len_unit.^5;
% Ke = 1/3*pi*density*a^3*u0^2*len_unit^5;
Ke = 7/23*Kv;
% K = 10/7*pi*density*a^3*u0^2*len_unit^5;
K = 30/23*Kv;