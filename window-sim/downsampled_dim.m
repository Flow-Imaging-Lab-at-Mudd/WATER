function Nd = downsampled_dim(dim, windows, overlap)
% Computes a downsampled dimension given the original dimension and
% windowing parameters.

% Compute the uniform dimension of downsampled fields.
Nd = 1 + floor((dim-windows)./(floor((1-overlap)*windows)));
disp('Dimensions of downsampled fields:')
disp(Nd)