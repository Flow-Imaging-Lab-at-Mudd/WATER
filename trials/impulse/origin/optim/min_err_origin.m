function error = min_err_origin(origin, smoother, num_ite, vf, props, I0)
% Wrapper objective function for the noisy impulse error with the choice
% of origin as the decision variable. The error is the average over a
% number of trials, as specified.
%
% Derek Li, December 2021

error = 0;
for i = 1: num_ite
    switch smoother
        case 'none'
            [err, ~, ~, ~, ~, ~, ...
                ~, ~, ~] = impulse_err_stats(vf, props, origin, I0);
        case 'box'
            [~, ~, err, ~, ~, ~, ...
                ~, ~, ~] = impulse_err_stats(vf, props, origin, I0);
        case 'gaussian'
            [~, ~, ~, ~, err, ~, ...
                ~, ~, ~] = impulse_err_stats(vf, props, origin, I0);
    end
    % Magnitude of error.
    error = error + norm(err);
end
% Average error over trials.
error = error / num_ite;