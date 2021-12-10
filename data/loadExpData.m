function vfs = loadExpData(status, frames)
% Load experimental data of turbulent vortex ring as velocity field
% objects given the time frames desired. Note that the frame indices given
% have an implied increment of 4, for these first few frames in the folder
% are invalid.

if ~exist('frames', 'var')
    frames = 1;
end

if ~isvector(frames)
    error('1D array of frame indices expected!')
end

% Load data from folder.
fdr = strcat(rootFolder, '\data\L18_Run1\');
data = dir(fdr);

% Subfolder names corresponding to frames.
frame_fdrs = strings(1, size(data, 1));
% Import subfolders in sorted order.
for i = 1: size(data, 1)
    frame_fdrs(i) = strcat(data(i).name, '\');
end

% Sort and select desired frames.
frame_fdrs = sort(frame_fdrs);
% The 5th frame is the first valid velocity dataset.
% First two entries not valid subfolders. Next two with null velocity
% fields.
frame_fdrs = frame_fdrs(5:end);

% Array of velocity fields, corresonding to different time frames.
vfs = cell(1, length(frames));

for i = 1: length(frames)
    switch status
        case 'raw'
            % Raw data.
            load(strcat(fdr, frame_fdrs(frames(i)), '\3DPIV_results.mat'), ...
                'xw', 'yw', 'zw', 'uw', 'vw', 'ww')
            vfs{i} = VelocityField.import_grid_separate(xw,yw,zw,uw,vw,ww);
        case 'post-processed'
            % Post-processed data.
            load(strcat(fdr, frame_fdrs(frames(i)), '\3DPIV_postprocessed_results_calibrated.mat'), ...
                'x', 'y', 'z', 'u', 'v', 'w')
            vfs{i} = VelocityField.import_grid_separate(x,y,z,u,v,w);
        otherwise
            error('Unsupported experimental data requested!')
    end
end

% Convenience for a single frame.
if length(frames) == 1
    vfs = vfs{1};
end
