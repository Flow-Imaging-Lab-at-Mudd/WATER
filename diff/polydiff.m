function F = polydiff(n)
% Derive finite difference formulas by interpolating 'n' uniformly spaced
% points to obtain a n-1 order polynomial, which is then differentiated.
% F[n-1, n, n] consists of differentiation upon to the n-1 order,
% difference formula at each of the n pivots, and the correspoinding n
% coefficients of y/h^k, in ascending order.

n = 6;
% Order of polynomial.
m = n - 1;
% Powers of polynomial.
p = (0: m)';

syms h
syms y [n 1]

x = transpose(0: h: m*h);

% Rows of distinct x values. Columns of increasing powers.
V = fliplr(vander(x));
% Polynomial coefficients.
a = V \ y;
% Construct polynomial.
syms x
P(x) = transpose(a)*x.^p;

% Compute coefficients of finite difference formulae pivoted at different
% indices.
F = zeros(m, n, n);

% Order of differentiation.
for k = 1: m
    dkP(x) = diff(P(x), x, k);
    % Pivot index.
    for i = 1: n
        dkP_i = dkP((i-1)*h);
        % Variable index.
        for j = 1: n
            F(k, i, j) = diff(dkP_i, y(j))*h^k;
        end
    end
end
