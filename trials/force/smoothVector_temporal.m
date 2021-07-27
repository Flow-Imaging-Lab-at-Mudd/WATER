function [v_s, dv_dt, spls] = smoothVector_temporal(v, t)
% Smooth vector quantity component-wise over time using a cubic
% spline.

v_s = zeros(size(v));
dv_dt = zeros(size(v));

% Splines fitted.
spls = cell(1, 3);

for d = 1: 3
    % Apply cubic smoothing spline. First try automatic
    % selection. If fails, require manual selection.
    vd = v(d,:);
    try
        [spls{d}, v_s(d,:), dv_dt(d,:), ~, ~, ~, ~, ~] = smoothspline(t, vd, 'auto', 3);
    catch
        [spls{d}, v_s(d,:), dv_dt(d,:), ~, ~, ~, ~, ~, figs] = smoothspline(t, vd, 'manual', 3);
        % Close all generated figures.
        close(figs)
    end
end
end