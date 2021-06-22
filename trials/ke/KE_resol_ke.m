% Plot KE error of vortical region (without random sampling) as KE noise.
sps = 0.01: 0.05: 1;
sps_count = size(sps, 2);
fr = 1;

% Introduce noise proportionally.
props = 0: 0.1: 7;
props_count = size(props, 2);

img_fdr = 'C:\Users\derek\flow\trials\ke\resol\global-ke\';
mkdir(img_fdr);


for i = 1: sps_count
    sp = sps(i);
    % Construct Hill vortex with specified resolution.
    [x, y, z, u, v, w, Mag] = hill_vortex_3D(sp, fr, 1, 1);
    vf = VelocityField.import_grid_separate(x,y,z,u,v,w);
    
    % Subtract freestream velocity to focus on central feature region.
    vf.addVelocity(-vf.U(1,1,1,:))
    
    KE_err_run(vf, props, fr);
    title(strcat('Feature Resolution $\frac{r}{s} = $', {' '}, string(fr/sp)))
    saveas(gcf, strcat(img_fdr, 'ke-', string(sp), 's.jpg'))
    close
end