function plt = plotVF(X, Y, range)
% range is a 3 x 2 array where each row records the start and end (inclusive)
% indices on the grid that corresponds to the rectangular region to be
% plotted.
%

if ndims(X) ~= 4 || ndims(Y) ~= 4
    error('4D formatted vector field required!')
end

if ~exist('range', 'var')
    range = [1 size(X,1); 1 size(X,2); 1 size(X,3)];
end

% Subset region plotted.
X = X(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :);
Y = Y(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :);

plt = figure;
quiver3(X(:,:,:,1), X(:,:,:,2), X(:,:,:,3), Y(:,:,:,1), Y(:,:,:,2), Y(:,:,:,3))
