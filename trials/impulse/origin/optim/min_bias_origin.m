function bias = min_bias_origin(origin, smoother, vf, props, I0)
% Wrapper objective function for the baseline impulse error with the choice
% of origin as the decision variable. The 'props' input is ignored for
% computing baseline error.
%
% Derek Li, December 2021

% Which situation is considered.
switch smoother
    case 'none'
        [~, ~, ~, ~, ~, ~, ...
            ~, ~, bias] = impulse_err_stats(vf, props, origin, I0);
    case 'box'
        [~, ~, ~, ~, ~, ~, ...
            bias, ~, ~] = impulse_err_stats(vf, props, origin, I0);
    case 'gaussian'
        [~, ~, ~, ~, ~, ~, ...
            ~, bias, ~] = impulse_err_stats(vf, props, origin, I0);
end
% Magnitude of error.
bias = norm(bias);