[x, y, z, u, v, w] = Hill_Vortex(0.25, 1, 1, 1, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

vf.clearNoise();
[~,~,sdu] = vf.noise_uniform_ugscaled(1);
% Flatten the matrices to compute correlation.
sdu = reshape(sdu, [], 1);
% Covariance matrix.
S = zeros(length(sdu), length(sdu));
for i = 1: length(sdu)
    for j = 1: i
        s = sdu(i)*sdu(j)*vf.corcoeff_tri(16, 0.75, i, j);
        if ~isreal(s)
            disp('unreal!')
            disp(s)
        end
        S(i,j) = s;
        if i~=j
            S(j,i) = s;
        end
    end
end
% L = chol(S);
% % Generate noise.
% N = zeros(size(vf.X_e));