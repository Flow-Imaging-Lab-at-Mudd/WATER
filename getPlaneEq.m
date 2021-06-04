function [orth, base, onPlane] = getPlaneEq(X)
% [orth, base, onPlane] = getPlaneEq(X)
%
% Given three positions stored as columns in the X matrix, compute a
% representation for a plane: the unit normal vector, the base position,
% and a boolean handle that decides whether a given point is on the plane.

% Tolerance for rounding error.
tol = 1e-4;

% Pick first position as base and use displacements from it to compute a
% unit normal vector.
orth = cross(X(:,2) - X(:,1), X(:,3) - X(:,1));
orth = orth / norm(orth);

base = X(:,1);

% Handle applicable exclusively to 4D matrix.
onPlane = @onPlane4d;
% (x) abs(dot(x - base, orth)) < tol;
    
    function on = onPlane4d(X)
        % Dimensionalize base and normal vectors for subtraction from 4D X.
        dim = size(X);
        dim = dim(1:3);
        Base = zeros(size(X));
        Base(:,:,:,1) = repmat(base(1), dim);
        Base(:,:,:,2) = repmat(base(2), dim);
        Base(:,:,:,3) = repmat(base(3), dim);
        Orth = zeros(size(X));
        Orth(:,:,:,1) = repmat(orth(1), dim);
        Orth(:,:,:,2) = repmat(orth(2), dim);
        Orth(:,:,:,3) = repmat(orth(3), dim);
        
        onplane = abs(dot(X - Base, Orth, 4)) < repmat(tol, dim);
        % Expand this logical 3D matrix indicating membership on plane to
        % 4D for direct point-wise multiplication with X.
        on = zeros(size(X));
        on(:,:,:,1) = onplane;
        on(:,:,:,2) = onplane;
        on(:,:,:,3) = onplane;
    end
end
