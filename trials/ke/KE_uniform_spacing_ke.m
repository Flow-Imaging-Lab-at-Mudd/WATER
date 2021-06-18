% Spacings of velocity.
sps = 0.01: 0.05: 1;
fr = 1;

for sp = sps
    KE_uniform
    title(strcat('Normalized Spacing:', {' '}, string(sp)))
    saveas(gcf, strcat('C:\Users\derek\flow\trials\ke\resol\global-ke\mild-noise\unit\0.66-1L\', 's=', string(sp), '.jpg'))
    close
end