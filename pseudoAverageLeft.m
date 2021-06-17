function V_pa = pseudoAverageLeft(V, winsize)
% V is presumed 3D matrix whose three dimensions are according to the
% meshgrid format. Only the left vertices of windows are valid in the
% matrix returned, which is of the same dimensions of the given.

% Flip x y dimension for meshgrid.
window = [winsize(2) winsize(1) winsize(3)];
% Pseudo-averaged matrix whose left indices alone are valid averages.
V_pa = zeros(size(V));

count = 0;
% Windowing average for the left indices.
for j = 0: window(1)-1
    for i = 0: window(2)-1
        for k = 0: window(3)-1
            V_pa(1: end-j, 1: end-i, 1: end-k) = ...
                V(1+j: end, 1+i: end, 1+k, end) + V_pa(1: end-j, 1: end-i, 1: end-k);
        end
    end
end

V_pa = V_pa / prod(window);
