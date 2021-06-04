function [orth, base, onPlane] = getPlaneEq(X)
% [orth, base, onPlane] = getPlaneEq(X)
%
% Given three positions stored as columns in the X matrix, compute a
% representation for a plane: the unit normal vector, the base position,
% and a boolean handle that decides whether a given point is on the plane.

% Pick first position as base and use displacements from it to compute a
% unit normal vector.
orth = cross(X(:,2) - X(:,1), X(:,3) - X(:,1));
orth = orth / norm(orth);

base = X(:,1);

% Rounding error for integral indices on grid handled automatically?
onPlane = @(x) dot(x - base, orth) == 0;
