% Add that folder plus all subfolders to the path.
global rootFolder
rootFolder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(rootFolder));
