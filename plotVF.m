function plt = plotVF(X, Y, scale, range)
% Customized plotting funnction for vector fields organized as 4D matrices.
% X is the matrix of positions on the grid, and Y the corresponding field.
% Note that the meshgrid indices of a 4D matrix is (y, x, z).
%
% range is a 3 x 2 array where each row records the start and end (inclusive)
% indices on the grid that corresponds to the rectangular region to be
% plotted.
%
% Note this plotting zooms in on the rectangular region specified by range.

if ndims(X) ~= 4 || ndims(Y) ~= 4 || ~isequal(size(X), size(Y))
    error('4D formatted vector field required!')
end

if ~exist('range', 'var')
    range = [1 size(X,1); 1 size(X,2); 1 size(X,3)];
end

% Subset region plotted.
if ~isequal(size(X), range * [-1 1]' + 1)
    X = X(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :);
    Y = Y(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :);
end

plt = figure;
quiver3(X(:,:,:,1), X(:,:,:,2), X(:,:,:,3), Y(:,:,:,1), Y(:,:,:,2), Y(:,:,:,3), scale)
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
