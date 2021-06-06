function on = skewPlaneMatrix(X, orth, base)
% SKEWPLANEMATRIX(X, orth, base) returns a 4D matrix corresponding to the
% positions in X, indicating whether positions therein are on the plane
% defined by the orthogonal vector and the base position given.

% Tolerance for rounding error.
tol = 1e-4;

% Upper-dimensionalize base and normal vectors for subtraction from 4D X.
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