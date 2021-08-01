% Generate distribution of unfiltered impulse error.

[x, y, z, u, v, w, Mag] = hill_vortex_3D(0.05, 1, 1, 1);
vf = VelocityField.import_grid_separate(x,y,z,u,v,w);

num_ite = 1;
err = zeros(3, num_ite);

origin = [0.21 -0.9 1.42]';


u_mean = vf.meanSpeed(0, 0);

% Normalizing factor for integrand.
i0 = mean(sqrt(sum(cross(VelocityField.sum3Vector(vf.X_e, -origin), vf.vort_e, 4).^2, 4)), 'all');

dims = [1 2 3];
% Iterate for a number of times.
for i = 1: num_ite
    % Noise profile.
    du = vf.noise_uniform(u_mean);
    % Compute integrand of impulse error without units scaling.
    dimp = cross(VelocityField.sum3Vector(vf.X_e, -origin), vf.curl(du), 4);
    for d = dims
        figure;
        histogram(du(:,:,:,d) / u_mean, 'Normalization', 'probability');
        xlabel(strcat('$\frac{\delta u_', vf.dim_str(d), '}{\bar{u}}$'))
        title(strcat('Distribution of $\delta u_', vf.dim_str(d), '$'))
        
        figure;
        histogram(dimp(:,:,:,d) / i0, 'Normalization', 'probability');
        xlabel(strcat('$\frac{|\delta i_', vf.dim_str(d), '|}{\bar{i}}$'))
        title(strcat('Distribution of $\delta i_', vf.dim_str(d), '$'))
    end
    
    % Distribution of noise and error magnitudes.
    figure;
    mag_du = sqrt(sum(du.^2, 4));
    histogram(mag_du / u_mean, 'Normalization', 'probability');
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    title(strcat('Distribution of $\delta u$'))
    
    figure;
    mag_dimp = sqrt(sum(dimp.^2, 4));
    histogram(mag_dimp / i0, 'Normalization', 'probability');
    xlabel('$\frac{|\delta i|}{\bar{i}}$')
    title(strcat('Distribution of $\delta i$'))
    
    pause
end
