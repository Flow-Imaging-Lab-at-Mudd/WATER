% The error trials for impulse are organized as follows, whose structrure
% is analogous to those of KE for noise propagation and resolution. By
% default, all tests are performed on the synthetic Hill's vortex.
%
% 'impulse_err_run.m' and 'impulse_err_stats.m' perform the basic
% operations of adding uniform random noise at different levels to the
% system, the former producing varieties of graphs and computing errors,
% the latter, calling the former, returning only statistics of errors,
% means and deviations. If a list of noise levels is feeded to the latter,
% an average will be computed, so we may prescribe a single level of noise
% and then iterate over a number of calls to 'impulse_err_stats' to account
% for stochasticity.
%
% The resolution test is performed in 'impulse_resol.m', which varies the
% resolution in the construction of Hill's vortex while holding the
% vortical radius constant. At each resolution, error statistics are
% computed from 'impulse_err_stats'. Again, we may desire a constant level
% of noise to generate an error vs. resolution graph. The resolution
% specified as parameter to 'impulse_error' to is feature resolution, r/s,
% how many vectors we have per feature, which is computed as feature length
% over global spacing. A minimum of 1 is required.
%
% To demonstrate that feature resolution, instead of say global resolution,
% is the determining factor to error, 'impulse_feature_resol.m' varies the
% feature radii while producing the same feature resolutions to generate
% a plot showing variations in the radii employed.
%
% Derek Li, September 2021
