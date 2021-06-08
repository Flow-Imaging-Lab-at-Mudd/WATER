function on = skewPlaneMatrix(X, orth, base, sz)
% SKEWPLANEMATRIX(X, orth, base) returns a matrix corresponding to the
% positions in X, indicating whether positions therein are on the plane
% defined by the orthogonal vector and the base position given.
%
% sz is the size of the last dimension used for element-wise multiplication.
% For velocity, sz = 4; for a scalar field, say error, sz = 1.

% Tolerance for rounding error.
tol = 1e-4;

% Upper-dimensionalize base and normal vectors for subtraction from 4D X.
dims = size(X);
dims = dims(1:3);

Base = dimen(base, dims);
Orth = dimen(orth, dims);

onplane = abs(dot(X - Base, Orth, 4)) < repmat(tol, dims);
% Expand this logical 3D matrix indicating membership on plane to
% the size of the last dimension for direct point-wise multiplication.
on = zeros([dims sz]);
for i = 1: sz
    on(:,:,:,i) = onplane;
end