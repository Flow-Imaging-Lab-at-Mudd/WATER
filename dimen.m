function V = dimen(v, dims)
% Replicate vector v to have preceding dimensions as specified by dims(1) x
% dims(2) x dims(3) x length(v). Size of dimension is upposed to be 3.

V = zeros([dims length(v)]);

for i = 1: length(v)
    V(:,:,:,i) = repmat(v(i), dims);
end
