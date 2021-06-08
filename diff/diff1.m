function nder = diff1(V, ind_inc, step, mode)
% Computes finite rate of change between two indices.

shifted = zeros(size(V));

ind_step = sum(ind_inc);

switch mode
    case 'right'
        switch find(ind_inc)
            % Shift according to direction and increment in index.
            case 1
                shifted(:, end-ind_step+1: end, :, :) = NaN;
                shifted(:, 1: end-ind_step, :, :) = V(:, 1+ind_step: end, :, :);
            case 2
                shifted(end-ind_step+1: end, :, :, :) = NaN;
                shifted(1: end-ind_step, :, :, :) = V(1+ind_step: end, :, :, :);
            case 3
                shifted(:, :, end-ind_step+1: end, :) = NaN;
                shifted(:, :, 1: end-ind_step, :) = V(:, :, 1+ind_step: end, :);
        end
    case 'left'
        switch find(ind_inc)
            case 1
        end
end

nder = (shifted - V) / step;