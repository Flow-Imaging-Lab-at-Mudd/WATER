function [dI, dId, dI_box, dI_gss, bias_box, bias_gss, vfd] = ...
    impulse_err_window_run(vf0, props, origin, I0, winsize, overlap)
% Incorporate the loss of resolution due to downsampling to the error 
% calculation. The accurate velocity field is first downsampled, and then
% noise is introduced to the downsampled field. The error is the difference
% between the impulse computed from the downsampled velocity after
% smoothing with the theoretical value of the original.
%
% 'I0' is the given theoretical or expected momentum taken as true.
%
% Derek Li, November 2021

% Theoretical momentum.
i0 = norm(I0);

% Create downsampled field and compute pure downsampling error.
vfd = vf0.downsample(winsize, overlap, 0);
% Impulse error from downsampling.
Id = vfd.impulse(origin, 0);
dId = (Id - I0) / i0;

% Introduce error to the downsampled field. The positions are maintained
% during downsampling. Thus the theoretical impulse used as reference level
% in impulse_err_run.m is accurate.
[dI, dI_box, dI_gss, bias_box, bias_gss] = impulse_err_run(vfd, props, origin, I0, false);