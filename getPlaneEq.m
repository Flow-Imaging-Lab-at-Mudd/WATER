function [orth, base, onPlane] = getPlaneEq(x)
% [orth, base, onPlane] = getPlaneEq(X)
%
% Given three positions stored as columns in the x matrix, compute a
% representation for a plane: the unit normal vector, the base position,
% and a boolean handle that decides whether a given point is on the plane.

% Pick first position as base and use displacements from it to compute a
% unit normal vector.
orth = cross(x(:,2) - x(:,1), x(:,3) - x(:,1));
orth = orth / norm(orth);

base = x(:,1);

% Handle applicable exclusively to 4D matrix.
onPlane = @(X) skewPlaneMatrix(X, orth, base);
