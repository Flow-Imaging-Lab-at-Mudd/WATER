function D = centricMatrix(num_points, diff_order, err_order)
% Construct a nth order differentiation matrix for a set of uniformly
% spaced points. Beware that, for generality, the differentiation matrix is
% not scaled by the proper power of step size, which should be multiplied
% once this differentiation matrix is obtained.

% Number of points needed to interpolate a polynomial.
n = diff_order + err_order;
if n > num_points
    error(strcat('Insufficient number of grid points to construct', ... 
        'differentiation matrix of the specified orders.', ' ', string(n), 'required!'))
end

% Compute finite difference formulas.
F = polydiff(n);
% Extract matrix for desired order of differentiation.
F = squeeze(F(diff_order, :, :));

% Range of indices where central formula applies. Note that for even
% cases, the "central" formula corresponds to the point n/2 + 1.
idx_i = floor(n/2) + 1;
% For even 'n', the last index is one farther right.
idx_f = num_points - (idx_i - 1) + (mod(n,2)==0)*1;
center_count = idx_f - idx_i + 1;

% Construct diagonal section of central coefficients.
C = spdiags(repmat(F(idx_i,:), center_count, 1), 0: n-1, center_count, num_points);
% Append coefficients for left and right boundary points.
L = [F(1:idx_i-1, :) zeros(idx_i-1, num_points-n)];
R = [zeros(num_points-idx_f, num_points-n) F(idx_i+1:end, :)];

D = [L; C; R];
