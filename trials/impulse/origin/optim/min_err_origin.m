% Objective function for minimizing the mean error (bias) by the
% selection of origin.

function err = min_err_origin(origin, smoother, vf, props, fr, u0)
    % Which situation is considered.
    switch smoother
        case 'none'
            [err, ~, ~, ~, ~, ~, ...
                ~, ~, ~] = impulse_err_stats(vf, props, origin, fr, u0);
        case 'box'
            [~, ~, err, ~, ~, ~, ...
                ~, ~, ~] = impulse_err_stats(vf, props, origin, fr, u0);
        case 'gaussian'
            [~, ~, ~, ~, err, ~, ...
                ~, ~, ~] = impulse_err_stats(vf, props, origin, fr, u0);
    end
    % Magnitude of error.
    err = norm(err);
    err
    origin
end