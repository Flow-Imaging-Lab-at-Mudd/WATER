function y = polyplot(p, x, col, sign)
% Plot a polynomial returned as vector by the polyfit function.
% Returns the values of the polynomial evaluated at the given points.

if ~exist('sign', 'var')
    sign = 1;
end

y = 0;
deg = length(p) - 1;

for k = 0: deg
    y = y + p(k+1) * x.^(sign*(deg-k));
end

plot(x, y, col)
% str = sprintf('$%.3f \epsilon_u^2 + (%f.3) \epsilon_u + (%f.3)$', p(1), p(2), p(3));