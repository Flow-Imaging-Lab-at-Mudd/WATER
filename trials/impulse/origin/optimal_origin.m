function [origin, err] = optimal_origin(vf, I0, origin0, ps_range)
% Shorthand for finding an origin that yields minimum error in impulse
% computation with respect to the given theoretical value 'I0' on the
% velocity field 'vf'. Noise is assumed to be added to 'vf' already in its
% field vf.N. The minimization algorithm used is 'patternsearch', which
% considers problems with stochasticity.
% 
% 'origin0' is the initial guess, and 'ps_range' specifies the range of
% origins to be considered, as a 3 x 2 matrix specifying the lower and
% upper bounds.

% Compute optimal origin for this specific noise profile.
[origin, err] = patternsearch(@(o) norm(vf.impulse(o, 1)-I0)/norm(I0), origin0, ...
    [],[],[],[], ps_range(:,1), ps_range(:,2));
