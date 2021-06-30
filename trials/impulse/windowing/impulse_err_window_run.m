function [dI, dId, dI_box, dI_gss, bias_box, bias_gss, vfd] = ...
    impulse_err_window_run(vf0, props, origin, fr, u0, winsize, overlap)
% Incorporate the loss of resolution due to downsampling to the error 
% calculation. The accurate velocity field is first downsampled, and then
% noise is introduced to the downsampled field. The error is the difference
% between the impulse computed from the downsampled velocity after
% smoothing with the theoretical value of the original.
%
% Derek Li, June 2021

% Theoretical momentum.
I0 = vf0.fluid.density*[0 2*pi*fr^3*u0*vf0.scale.len^4 0]';
i0 = norm(I0);

% Create downsampled field and compute pure downsampling error.
vfd = vf0.downsample(winsize, overlap, 0);
% Shift between noise-free KE vals.
Id = vfd.impulse(0, origin);
% KE error from downsampling.
dId = (Id - I0) / i0;

% Introduce error to the downsampled field. The positions are maintained
% during downsampling. Thus the theoretical impulse used as reference level
% in impulse_err_run.m is accurate.
[dI, dI_box, dI_gss, bias_box, bias_gss] = impulse_err_run(vfd, props, origin, fr, u0);