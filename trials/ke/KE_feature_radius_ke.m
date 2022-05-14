% Feature radius.
radii = 0.1: 0.1: 1;
radii_count = length(radii);
sp = 0.1;

% Proportions of noise to mean speed.
props = 0: 0.1: 3;

% Containers for data over iterations.
dK = zeros(length(props), radii_count);
dK_box = zeros(length(props), radii_count);
dK_gss = zeros(length(props), radii_count);
bias_box = zeros(radii_count);
bias_gss = zeros(radii_count);

img_fdr = strcat('C:\Users\derek\flow\trials\ke\feature\global-ke\s=', string(sp), '\');
mkdir(img_fdr);

for i = 1: radii_count
    fr = radii(i);
    [x, y, z, u, v, w, ~] = hill_vortex_3D(sp, fr, 1, 1);
    vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
    % Focus on the global region with the freestream velocity subtracted.
    vf.addVelocity(-vf.U(1,1,1,:));
    % Run script for KE error sampling.
    [dK(:,i), dK_box(:,i), dK_gss(:,i), bias_box(i), bias_gss(i)] = ...
        KE_err_run(vf, props, fr);
    
    K_noise = dK(:, i);
    
    figure;
    scatter(K_noise, abs(dK_box(:, i)), 'r', 'filled')
    hold on
    scatter(K_noise, abs(dK_gss(:, i)), 'b', 'filled')
    hold on

    % 1-1 line.
    plot(K_noise, K_noise, 'black')
    hold on

    % Smoother biases.
    mean_bias_box = mean(abs(bias_box(i)));
    sd_bias_box = std(abs(bias_box(i)));
    yline(mean_bias_box, '-', 'Color', 'r')

    hold on
    mean_bias_gss = mean(abs(bias_gss(i)));
    sd_bias_gss = std(abs(bias_gss(i)));
    yline(mean_bias_gss, '-', 'Color', 'b')

    legend('box-filtered', 'Gaussian-filtered', ...
        '$y=x$ identity line', ...
        strcat('box bias $\kappa = $', string(mean_bias_box), '$\pm$', string(sd_bias_box)), ...
        strcat('Gaussian bias $\kappa = $', string(mean_bias_gss), '$\pm$', string(sd_bias_gss)), ...
        'Interpreter', 'latex')
    xlabel('Unfiltered $\left|\frac{\delta K}{K}\right|$')
    ylabel('Filtered $\left|\frac{\delta K}{K}\right|$')
    title(strcat('Feature Resolution $\frac{r}{s} = $', {' '}, string(fr/sp)))
    saveas(gcf, strcat(img_fdr, 'ke-', string(fr), 'r', '.jpg'))
    close
end