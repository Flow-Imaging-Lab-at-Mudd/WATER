% This function simulates 3DPIV windowing effects (i.e. converts from voxel
% to interrogation window scale)
% Written by Leah Mendelson & Derek Li

function [xwin,ywin,zwin,uwin,vwin,wwin] = PIV_window_sim(x,y,z,u,v,w,winsize,overlap)

    % inputs
    % x,y,z: input vector field coordinates (analogous to voxels)
    % u,v,w: velocity vector components at input field resolution
    % winsize: simulated interrogation window size, specified as a power of
    % 2 (typical range 16-128 voxels). Code assumes a uniform interrogation
    % window size in each direction
    % overlap: Overlap between interrogation windows as a fraction of the window
    % size in voxels (typically 0.5 or 0.75). Code assumes a uniform
    % overlap in each direction
    %
    % outputs
    % xwin,ywin,zwin: output vector field coordinates (in pixels, analogous
    % to interrogation window resolution)
    % uwin,vwin,wvin: filtered velocity vector components at output field
    % resolution
    
    if isequal(size(winsize), [1  1])
        winsize = winsize*ones(1, 3);
    end
    
    % Dimensions of grid.
    xdim = size(x, 2);
    ydim = size(x, 1);
    zdim = size(x, 3);
    
    % Box averaging of velociy values.
    % after this step, vectors are still on input grid, but filtered like
    % PIV processing; only those on left vertices are valid.
    ufilt = pseudoAverageLeft(u, winsize);
    vfilt = pseudoAverageLeft(v, winsize);
    wfilt = pseudoAverageLeft(w, winsize);
        
    % generate coordinates for each interrogation window (center of window)
    if min(winsize) >= 1
        % Frame shifts of indices.
        shift_x = winsize(1)*(1 - overlap(1));
        shift_y = winsize(2)*(1 - overlap(2));
        shift_z = winsize(3)*(1 - overlap(3));
        
        % Dimensions of downsampled grid.
        xdim_d = floor((xdim - winxsize(1)) / shift_x) + 1;
        ydim_d = floor((ydim - winxsize(2)) / shift_y) + 1;
        zdim_d = floor((zdim - winxsize(3)) / shift_z) + 1;
        
        % edge handling
        % last coordinate for which enough data is available for a complete interrogation window of size boxes(dc)
        endIndex_x = shift_x * (xdim_d-1) + 1;
        endIndex_y = shift_y * (ydim_d-1) + 1;
        endIndex_z = shift_z * (zdim_d-1) + 1;
        
        % Averaging.
        
        % Downsample. 
        xrange = 1: shift_x: endIndex_x;
        yrange = 1: shift_y: endIndex_y;
        zrange = 1: shift_z: endIndex_z;
        
        uwin = ufilt(xrange, yrange, zrange);
        vwin = vfilt(xrange, yrange ,zrange);
        wwin = wfilt(xrange, yrange ,zrange);
        
        xwin = ufilt(xrange, yrange, zrange);
        ywin = vfilt(xrange, yrange ,zrange);
        zwin = wfilt(xrange, yrange ,zrange);
        
        [xwin, ywin, zwin] = meshgrid(xx,yy,zz);
    else
        disp('Error. Window size must be larger than 1 voxel')
        return
    end
        
end