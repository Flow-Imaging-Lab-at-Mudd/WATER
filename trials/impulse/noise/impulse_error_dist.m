% Generate distribution of unfiltered error quantities in impulse
% computation.
%
% Derek Li, August 2021

[x, y, z, u, v, w, Mag] = Hill_Vortex(0.05, 1, 1, 1, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

% Computed impulse error to generate distribution.
num_ite = 10000;
err = zeros(3, num_ite);

origin = [0.21 -0.9 1.42]';

u_mean = vf.meanSpeed(0, 0);

% Integrand of impulse.
inte = cross(VelocityField.sum3Vector(vf.X_e, -origin), vf.vort_e, 4);
% Normalizing factor for integrand.
i0 = mean(sqrt(sum(inte.^2, 4)), 'all');

% Magnitude of vorticity for normalization.
vort_mag = sqrt(sum(vf.vort_e.^2, 4));
vort_mean = mean(vort_mag, 'all');

dims = [];
% Iterate for a number of times.
for i = 1: num_ite
    % Noise profile.
    du = vf.noise_uniform(u_mean);
    % Vorticity error.
    dvort = vf.curl(du);
    % Compute integrand of impulse error without units scaling.
    dimp = cross(VelocityField.sum3Vector(vf.X_e, -origin), vf.curl(du), 4);
    for d = dims
        % Dimensional values, normalized.
        dud = reshape(du(:,:,:,d), [], 1) / u_mean;
        dvortd = reshape(dvort(:,:,:,d), [], 1) / vort_mean;
        dimpd = reshape(dimp(:,:,:,d), [], 1) / i0;
        
        % Fit normal distributions.
        nd_dvortd = fitdist(dvortd, 'Normal');
        nd_dimpd = fitdist(dimpd, 'Normal');
        
        % Velocity noise distribution.
        figure;
        histogram(dud, 'Normalization', 'pdf');
        xlabel(strcat('$\frac{\delta u_', vf.dim_str(d), '}{\bar{u}}$'))
        title(strcat('Distribution of $\delta u_', vf.dim_str(d), '$'))
        
        % Vorticity error distribution.
        figure;
        histogram(dvortd, 'Normalization', 'pdf');
        hold on
        minval = min(dvortd);
        maxval = max(dvortd);
        range = minval: (maxval-minval)/100: maxval;
        plot(range, pdf(nd_dvortd, range), 'LineWidth', 2)
        xlabel(strcat('$\frac{|\delta \omega_', vf.dim_str(d), '|}{\bar{\omega}}$'))
        title(strcat('Distribution of $\delta \omega_', vf.dim_str(d), '$'))
        
        % Distribution of integrand.
        figure;
        histogram(dimpd, 'Normalization', 'pdf');
        hold on
        minval = min(dimpd);
        maxval = max(dimpd);
        range = minval: (maxval-minval)/100: maxval;
        plot(range, pdf(nd_dimpd, range), 'LineWidth', 2)
        xlabel(strcat('$\frac{|\delta i_', vf.dim_str(d), '|}{\bar{i}}$'))
        title(strcat('Distribution of $\delta i_', vf.dim_str(d), '$'))
    end
    
    % Distribution of velocity noise.
    figure;
    du_mag = reshape(sqrt(sum(du.^2, 4)), [], 1);
    histogram(du_mag / u_mean, 'Normalization', 'pdf');
    hold on
    % Fit normal distribution.
    nd_du = fitdist(du_mag, 'Normal');
    minval = min(du_mag);
    maxval = max(du_mag);
    range = minval: (maxval-minval)/100: maxval;
    plot(range, pdf(nd_du, range), 'LineWidth', 2)
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    title(strcat('Distribution of $\delta u$'))
    
    % Distribution of vorticity error magnitude.
    figure;
    dvort_mag = sqrt(sum(dvort.^2, 4));
    histogram(dvort_mag / vort_mean, 'Normalization', 'pdf');
    xlabel('$\frac{|\delta \omega|}{\bar{\omega}}$')
    title(strcat('Distribution of $\delta \omega$'))
    
    % Distribution of integrand error magnitude.
    figure;
    dimp_mag = sqrt(sum(dimp.^2, 4));
    histogram(dimp_mag / i0, 'Normalization', 'pdf');
    xlabel('$\frac{|\delta i|}{\bar{i}}$')
    title(strcat('Distribution of $\delta i$'))
    
    % Sum integrad to obtain impulse error.
    err(:, i) = sum(dimp, [1 2 3]);
    
    pause
end

% Magnitude of impulse used for normalizing. Unscale the units.
I_mag = norm(squeeze(sum(inte, [1 2 3])));

% Render impulse error distribution.
dims = [1 2 3];
for d = dims
    figure;
    errd = reshape(err(d, :)/I_mag, [], 1);
    histogram(errd, 'Normalization', 'pdf')
    hold on
    % Fit normal distribution.
    nd_errd = fitdist(errd, 'Normal');
    minval = min(errd);
    maxval = max(errd);
    range = minval: (maxval-minval)/100: maxval;
    plot(range, pdf(nd_errd, range), 'LineWidth', 2)
    xlabel(strcat('$\frac{\delta I_', vf.dim_str(d), '}{\bar{I}}$'))
    title('Distribution of impulse error')
end
