% Objective function for minimizing the baseline error (bias) by the
% selection of origin.

function bias = min_bias_origin(origin, smoother, vf, props, fr, u0)
    % Which situation is considered.
    switch smoother
        case 'none'
            [~, ~, ~, ~, ~, ~, ...
                ~, ~, bias] = impulse_err_stats(vf, props, origin, fr, u0);
        case 'box'
            [~, ~, ~, ~, ~, ~, ...
                bias, ~, ~] = impulse_err_stats(vf, props, origin, fr, u0);
        case 'gaussian'
            [~, ~, ~, ~, ~, ~, ...
                ~, bias, ~] = impulse_err_stats(vf, props, origin, fr, u0);
    end
    % Magnitude of error.
    bias = norm(bias);
end