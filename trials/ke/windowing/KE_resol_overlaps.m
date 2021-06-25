% Identical to KE_resol.m while plotting for different overlap poportions.

ops = [0.25 0.5 0.75];
ops_count = size(ops, 2);

% Desired error level.
err_level = 0.1;

% Constant feature radius used while varying global resolution.
fr = 1;

% Downsampling parameters.
winsize = 16;

% Maximum freature size tried.
max_fres = 20;
% Increment.
fres_inc = 1;

% Generate range of downsampled spacings for evenly spaced feature
% resolutions.
% Global spacing of downsampled data.
dsp = zeros(1, floor((max_fres-1)/fres_inc) + 1);
% Minimal feature resolution of 1.
dsp(1) = 1 / fr;

for i = 2: size(dsp, 2)
    dsp(i) = fr / (fres_inc + fr/dsp(i-1));
end

% Reverse ordering of spacing so that high resolutions precedes low.
dsp = flip(dsp);

% Compute initial resolutions used. Now this depends on the proportion of
% overlap.
resol = dsp' ./ min((1-ops)*winsize, winsize-1);
resol_count = size(resol, 1);

% Feature resolution.
fres = fr ./ dsp;

% Containers for data across all runs.
% Errors here are mean absolute errors.
err_box = NaN(resol_count, ops_count);
err_gss = NaN(resol_count, ops_count);

bias_box = NaN(resol_count, ops_count);
bias_gss = NaN(resol_count, ops_count);

% Downsampling bias.
dKd = NaN(resol_count, ops_count);

% Introduce noise proportionally.
props = 0: 0.1: 3;

for o = 1: ops_count
    overlap = ops(o);
    for k = 1: resol_count
        % Construct Hill vortex with specified resolution.
        sp = resol(k, o);
        [x, y, z, u, v, w, ~] = hill_vortex_3D(sp, fr, 1, 1);
        vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
        % Subtract freestream velocity to focus on central feature region.
        vf.addVelocity(-vf.U(1,1,1,:))
        
        try
        % Incorporate windowing.
        [~, dKd(k,o),  ~, ~, bias_box(k,o), bias_gss(k,o), ...
            err_box(k,o), err_gss(k,o)] = ...
            KE_err_window_run(vf, winsize, overlap, props);
        catch
            % Lower resolutions impossible.
            break
        end
    end
end

% Absolute errors.
smoother_bias_box = abs(bias_box);
smoother_bias_gss = abs(bias_gss);

fres_min_err = zeros(2, ops_count);
fres_min_bias = zeros(2, ops_count);


abscissa = 1 ./ fres;

% Superimpose plots for different overlap ratios.
for o = 1: ops_count
    % Lowest feature resolution to achieve the desired error level.
    try
        fres_min_err(1, o) = fres(end - find(flip(err_box(:,o)) < err_level, 1) + 1);
    catch
        fres_min_err(1, o) = -1;
    end
    try
        fres_min_err(2, o) = fres(end - find(flip(err_gss(:,o)) < err_level, 1) + 1);
    catch
        fres_min_err(2, o) = -1;
    end

    % Lowest feature resolution to achieve the desired bias level.
    try
        fres_min_bias(1, o) = fres(end - find(flip(smoother_bias_box(:,o)) < err_level, 1) + 1);
    catch
        fres_min_bias(1, o) = -1;
    end
    try
        fres_min_bias(2, o) = fres(end - find(flip(smoother_bias_gss(:,o)) < err_level, 1) + 1);
    catch
        fres_min_bias(2, o) = -1;
    end
end

% Colors in plots corresponding to 25%, 50%, and 75% errors.
cols = {'g', 'b', 'r'};

% Box filter mean errors.
figure;
% Desired error line.
yline(err_level, '-')
for o = 1: ops_count
    % Eliminate under-resolved data.
    err = err_box(:, o);
    if find(isnan(err), 1)
        idx = find(isnan(err), 1);
    else
        idx = resol_count + 1;
    end
    % Valid resolutions for this overlap.
    odsp = dsp(1:idx-1);
    ofres = fres(1:idx-1);
    oerr_box = err_box(1:idx-1, o);
    
    overlap = ops(o);
    
    hold on
    scatter(ofres, oerr_box, 'ko', 'MarkerFaceColor', cols{o}, 'LineWidth', 1)
end
legend({strcat(string(err_level*100), '\% error line')...
    '25\% overlap', '50\% overlap', '75\% overlap'}, ...    
        'Interpreter', 'latex')
xlabel(strcat('Downsampled feature resolution $\frac{r}{s}$'))
ylabel('$\bar{\left|\frac{\delta K}{K}\right|}$')
title(strcat('Mean Error of box Smoother over $\delta u = $', {' '}, ...
    string(props(1)), '-', string(props(end)*100), '\% $\bar{u}$'))


% Gaussian filter mean errors.
figure;
% Desired error line.
yline(err_level, '-')
for o = 1: ops_count
    % Eliminate under-resolved data.
    err = err_gss(:, o);
    if find(isnan(err), 1)
        idx = find(isnan(err), 1);
    else
        idx = resol_count + 1;
    end
    % Valid resolutions for this overlap.
    odsp = dsp(1:idx-1);
    ofres = fres(1:idx-1);
    obias_gss = bias_gss(1:idx-1, o);
    oerr_gss = err_gss(1:idx-1, o);
    
    overlap = ops(o);
    
    hold on
    scatter(ofres, oerr_gss, 'ko', 'MarkerFaceColor', cols{o}, 'LineWidth', 1)
end
legend({strcat(string(err_level*100), '\% error line')...
    '25\% overlap', '50\% overlap', '75\% overlap'}, ...    
        'Interpreter', 'latex')
xlabel(strcat('Downsampled feature resolution $\frac{r}{s}$'))
ylabel('$\bar{\left|\frac{\delta K}{K}\right|}$')
title(strcat('Mean Error of Gaussian Smoother over $\delta u = $', {' '}, ...
    string(props(1)), '-', string(props(end)*100), '\% $\bar{u}$'))


% Box filter bias.
figure;
% Desired error line.
yline(err_level, '-')
for o = 1: ops_count
    % Eliminate under-resolved data.
    err = err_box(:, o);
    if find(isnan(err), 1)
        idx = find(isnan(err), 1);
    else
        idx = resol_count + 1;
    end
    % Valid resolutions for this overlap.
    odsp = dsp(1:idx-1);
    ofres = fres(1:idx-1);
    obias_box = smoother_bias_box(1:idx-1, o);
    
    overlap = ops(o);
    
    hold on
    scatter(ofres, obias_box, 'ko', 'MarkerFaceColor', cols{o}, 'LineWidth', 1)
end
legend({strcat(string(err_level*100), '\% error line')...
    '25\% overlap', '50\% overlap', '75\% overlap'}, ...    
        'Interpreter', 'latex')
xlabel(strcat('Downsampled feature resolution $\frac{r}{s}$'))
ylabel('$\bar{\left|\frac{\delta K}{K}\right|}$')
title(strcat('Box Smoother Bias'))

% Gaussian filter bias.
figure;
% Desired error line.
yline(err_level, '-')
for o = 1: ops_count
    % Eliminate under-resolved data.
    err = err_box(:, o);
    if find(isnan(err), 1)
        idx = find(isnan(err), 1);
    else
        idx = resol_count + 1;
    end
    % Valid resolutions for this overlap.
    odsp = dsp(1:idx-1);
    ofres = fres(1:idx-1);
    obias_gss = smoother_bias_gss(1:idx-1, o);
    
    overlap = ops(o);
    
    hold on
    scatter(ofres, obias_gss, 'ko', 'MarkerFaceColor', cols{o}, 'LineWidth', 1)
end
legend({strcat(string(err_level*100), '\% error line')...
    '25\% overlap', '50\% overlap', '75\% overlap'}, ...    
        'Interpreter', 'latex')
xlabel(strcat('Downsampled feature resolution $\frac{r}{s}$'))
ylabel('$\bar{\left|\frac{\delta K}{K}\right|}$')
title(strcat('Gaussian Smoother Bias'))
