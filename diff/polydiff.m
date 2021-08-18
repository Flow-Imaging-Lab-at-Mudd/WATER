function [F, P] = polydiff(n)
% Derive symbolically finite difference formulas by interpolating 'n'
% uniformly spaced points to obtain a n-1 order polynomial, which is then
% differentiated. F[n-1, n, n] consists of differentiation upon to the n-1
% order, difference formula at each of the n pivots, and the correspoinding
% n coefficients of y/h^k, in ascending order.
%
% polydiff retains a table of pre-computed finite difference formulae for
% look-up, up to the order of 40. This table is stored in the global
% variable PDF.
%
% A symbolic handle of the interpolated polynomial P(x) is also returned.
% Since a symbolic function created in this scope does not automatically
% declare the symbolic variables it depends on, they must be declared
% manually when used. P is a function of 'x', with step size 'h' and values
% 'y1' ... 'yn-1' corresponding to 0 h 2h ... (n-1)h.
%
% For interpolating multidimensional ordinates, the derivative matrix is
% simply applied to each dimension, for the result is general.

% Order of polynomial.
m = n - 1;

% Maximum order of finite difference stored. This is only used for
% initialization and may be modified if a higher order formula is computed.
max_order = 40;

global PDF
% Check if the formula has already been computed.
if length(PDF) > 1
    % Expand the table if a higher order formula is demanded.
    if m > size(PDF, 1)
        PDF = [PDF; cell((floor(m/10)+1)*10, 2)];
    elseif ~isempty(PDF{m, 1})
        F = PDF{m, 1};
        P = PDF{m, 2};
        return
    end
else
    % Initialize table. First column 3D matrix, second column polynomial.
    PDF = cell(max_order, 2);
end

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
        % Ordinate index.
        for j = 1: n
            F(k, i, j) = diff(dkP_i, y(j))*h^k;
        end
    end
end

% % Store differentiation matrices and polynomial.
% PDF{m, 1} = F;
% PDF{m, 2} = P;
