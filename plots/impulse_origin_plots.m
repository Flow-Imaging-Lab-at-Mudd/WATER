% Load computed data of impulse errors at uniformly sampled impulse
% locations.
load('origin-graph.mat')

% Font.
font = 'Times New Roman';
fontSize = 10;
% Only for title.
fontWeight = 'normal';

t = tiledlayout(1, 2);

% Error magnitude after noise without smoothing.
nexttile;
vfp.plotScalar(mag_err(:,:,:,2), 0, '');
title('(a) Impulse error using different origins')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = '$\frac{|\delta\vec{I}|}{I}$';
c.Label.Rotation = 0;
c.Label.FontSize = 12;
% c.Limits = [0 0.08];

axis square

% Plot of standard deviations.
nexttile;
mag_err_sd = squeeze(mag_err_sd);
vfp.plotScalar(mag_err_sd(:,:,:,2), 0, '');
title('(b) Variation of impulse error over trials')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = 'SD$\left(\frac{|\delta\vec{I}|}{I}\right)$';
c.Label.Rotation = 0;
c.Label.FontSize = 12;
% c.Limits = [0 0.08];

axis square