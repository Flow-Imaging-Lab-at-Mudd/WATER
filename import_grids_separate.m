function [X, U] = import_grids_separate(xw, yw, zw, uw, vw, ww)

% Pack positions and velocities compactly as 3-vectors in extra dimension.
X = xw;
X(:,:,:,2) = yw;
X(:,:,:,3) = zw;

U = uw;
U(:,:,:,2) = vw;
U(:,:,:,3) = ww;
