% Repeat impulse_err_run.m for a number of iterations to account for
% stochastic effects using a Hill's vortex.

% Parameters of Hill's vortex.
r0 = 1;
sp = 0.05;
u0 = 1;

[x, y, z, u, v, w, Mag] = Hill_Vortex(sp, r0, u0, 1, 1);
vf = VelocityField.importCmps(x, y, z, u, v, w);

origin = [0 0 0]';
u_mean = vf.meanSpeed(0, 0);

% Proportions of error.
props = 0: 0.5: 3;
props_count = size(props, 2);

I0 = HillImpulse(vf.fluid.density, vf.scale.len, r0, u0);
i0 = I0(2);

% Number of iterates at each level of noise.
num_ite = 5;

% Error in impulse computation given noise.
dI = zeros(3, props_count, num_ite);
% Box smoothing.
dI_box = zeros(3, props_count, num_ite);
% Gaussian smoothing.
dI_gss = zeros(3, props_count, num_ite);

for k = 1: num_ite
    [dI(:,:,k), dI_box(:,:,k), dI_gss(:,:,k),~,~] = ...
        impulse_err_run(vf, props, origin, I0, 0);
end

di = squeeze(sqrt(sum(dI.^2, 1)));
di_box = squeeze(sqrt(sum(dI_box.^2, 1)));
di_gss = squeeze(sqrt(sum(dI_gss.^2, 1)));

abs_dI = abs(dI);
abs_dI_box = abs(dI_box);
abs_dI_gss = abs(dI_gss);

% Baseline smoother biases.
bias_box = dI_box(:, 1, 1);
bias_gss = dI_gss(:, 1, 1);

mag_bias_box = norm(bias_box);
mag_bias_gss = norm(bias_gss);

abs_bias_box = abs(bias_box);
abs_bias_gss = abs(bias_gss);

% Dimension, i.e., x, y, z, to plot, specified correspondingly by 1, 2, 3.
dims = [2];
dim_str = {'x', 'y', 'z'};

%%%%%%%%%%% Plot signed impulse error %%%%%%%%%%%%%

% Font size for plotting titles and axis labels.
fsize = 11;

for dim = dims
    % Unfiltered error.
    figure;
    for p = 1: props_count
        % Repeated noise level for plotting.
        pp = repmat(props(p), 1, num_ite);
        scatter(pp, squeeze(dI(dim, p, :)))
        hold on
    end
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    ylabel(strcat('$\frac{\delta I_', dim_str{dim}, '}{I}$'))
    title(strcat('Unfiltered $\hat{', dim_str{dim}, '}$ Impulse Error'))
    ax = gca;
    ax.FontSize = fsize;
    
    % Box-filtered.
    figure;
    for p = 1: props_count
        % Repeated noise level for plotting.
        pp = repmat(props(p), 1, num_ite);
        scatter(pp, squeeze(dI_box(dim, p, :)), 'r', 'filled')
        % Comparison with unfiltered.
        hold on
        scatter(pp, dI(dim, p, :))
        hold on
    end
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    ylabel(strcat('$\frac{\delta I_', dim_str{dim}, '}{I}$'))
    legend({'box-filtered', 'unfiltered'})
    title(strcat('Box-filtered $\hat{', dim_str{dim}, '}$ Impulse Error'))
    ax = gca;
    ax.FontSize = fsize;
    
    % Gaussian-filtered.
    figure;
    for p = 1: props_count
        % Repeated noise level for plotting.
        pp = repmat(props(p), 1, num_ite);
        scatter(pp, squeeze(dI_gss(dim, p, :)), 'b', 'filled')
        % Comparison with unfiltered.
        hold on
        scatter(pp, squeeze(dI(dim, p, :)))
        hold on
    end
    xlabel('$\frac{|\delta u|}{\bar{u}}$')
    ylabel(strcat('$\frac{\delta I_', dim_str{dim}, '}{I}$'))
    legend({'Gaussian-filtered', 'unfiltered'})
    title(strcat('Gaussian-filtered $\hat{', dim_str{dim}, '}$ Impulse Error'))
    ax = gca;
    ax.FontSize = fsize;
end

%%%%%%%%%%%% Magnitude Plots %%%%%%%%%%%%%
% Unfiltered error.
figure;
for p = 1: props_count
    % Repeated noise level for plotting.
    pp = repmat(props(p), 1, num_ite);
    scatter(pp, di(p, :))
    hold on
end
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('$\frac{|\delta I|}{\bar{I}}$')
title('Magnitude of unfiltered impulse error')
ax = gca;
ax.FontSize = fsize;

% Box-filtered.
figure;
for p = 1: props_count
    % Repeated noise level for plotting.
    pp = repmat(props(p), 1, num_ite);
    scatter(pp, di_box(p, :), 'r', 'filled')
    hold on
    scatter(pp, di(p, :))
    hold on
end
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('$\frac{|\delta I|}{\bar{I}}$')
legend({'box-filtered', 'unfiltered'})
title('Magnitude of box-filtered impulse error')
ax = gca;
ax.FontSize = fsize;

% Gaussian-filtered.
figure;
for p = 1: props_count
    % Repeated noise level for plotting.
    pp = repmat(props(p), 1, num_ite);
    scatter(pp, di_gss(p, :), 'b', 'filled')
    hold on
    scatter(pp, di(p, :))
    hold on
end
xlabel('$\frac{|\delta u|}{\bar{u}}$')
ylabel('$\frac{|\delta I|}{\bar{I}}$')
legend({'Gaussian-filtered', 'unfiltered'})
title('Magnitude of Gaussian-filtered impulse error')
ax = gca;
ax.FontSize = fsize;