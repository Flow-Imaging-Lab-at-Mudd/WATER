% This function simulates 3D PIV windowing effects (i.e. converts from voxel
% to interrogation window scale)
% Written by Leah Mendelson & Derek Li

function [Xwin, Uwin] = PIV_window_sim(X, U, winsize, overlap, Xscale)
% inputs
% X: 4D matrix of input vector field coordinates (analogous to voxels)
% U: Matching 4D matrix of velocity vector components at input field resolution
% winsize: simulated interrogation window size, specified as a power of
% 2 (typical range 16-128 voxels). Can be given a uniform scalar value
% or 3-vector of x y z window sizes.
% overlap: Overlap between interrogation windows as a fraction of the window
% size in voxels (typically 0.5 or 0.75). A uniform overlap in each direction
% is assumed.
% Xscale: Vector for how the subsetted positions are to be scaled, where
% all positions will be non-negative and proceed from (0 0 0) by the
% increments specified by 'Xscale'.
%
% outputs
% Xwin: output vector field coordinates (in pixels, analogous
% to interrogation window resolution)
% Uwin: filtered velocity vector components at output field
% resolution

if isequal(size(winsize), [1  1])
    winsize = winsize*ones(1, 3);
end

if isequal(size(overlap), [1  1])
    overlap = overlap*ones(1, 3);
end

% Dimensions of grid.
xdim = size(X, 2);
ydim = size(X, 1);
zdim = size(X, 3);

% Ensure overlap between windows is in integers.
overlap = round(winsize .* overlap);
overlap = min([overlap; winsize-1], [], 1);

% Box averaging of velociy values.
% after this step, vectors are still on input grid, but filtered like
% PIV processing.

ker = 1/prod(winsize) * ones(winsize);

Ufilt = convn(U, ker, 'same');

% Ufilt = pseudoAverageLeft(U, winsize);

% Properly subset the averaged velocities and create windowed position meshgrid.
if min(winsize) >= 1
    % Frame shifts of indices.
    shift_x = winsize(1) - overlap(1);
    shift_y = winsize(2) - overlap(2);
    shift_z = winsize(3) - overlap(3);

    % Dimensions of downsampled grid.
    xdim_d = floor((xdim - winsize(1)) / shift_x) + 1;
    ydim_d = floor((ydim - winsize(2)) / shift_y) + 1;
    zdim_d = floor((zdim - winsize(3)) / shift_z) + 1;

    % Edge handling
    % Last left coordinate for which enough data is available for a
    % complete interrogation window。
    endIndex_x = shift_x * (xdim_d-1) + 1;
    endIndex_y = shift_y * (ydim_d-1) + 1;
    endIndex_z = shift_z * (zdim_d-1) + 1;

    % Downsample. Left indices.
    xrange = 1: shift_x: endIndex_x;
    yrange = 1: shift_y: endIndex_y;
    zrange = 1: shift_z: endIndex_z;
    
    % Use central indices.
    Uwin = Ufilt(yrange + floor(winsize(2)/2 - 1), ...
        xrange + floor(winsize(1)/2 - 1), ...
        zrange + floor(winsize(3)/2 - 1), :);

    % If original coordinates are to be retained.
    if Xscale == 1
        [Xwin(:,:,:,1), Xwin(:,:,:,2), Xwin(:,:,:,3)] = ...
            meshgrid(X(1,xrange,1,1) + (X(1,winsize(1),1,1)-X(1,1,1,1))/2, ...
                X(yrange,1,1,2) + (X(winsize(2),1,1,2)-X(1,1,1,2))/2, ...
                X(1,1,zrange,3) + (X(1,1,winsize(3),3)-X(1,1,1,3))/2);
    else
        [Xwin(:,:,:,1), Xwin(:,:,:,2), Xwin(:,:,:,3)] = ...
            meshgrid((xrange-1)*Xscale(1), (yrange-1)*Xscale(2), (zrange-1)*Xscale(3));
    end
else
    error('Window size must be larger than 1 voxel!')
end
