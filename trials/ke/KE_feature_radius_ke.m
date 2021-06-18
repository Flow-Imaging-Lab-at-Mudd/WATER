% Feature radius.
frs = 0.1: 0.1: 1;
sp = 0.05;

for fr = frs
    KE_focus_feature
    title(strcat('Feature-Spacing Product $rs = $', {' '}, string(fr*sp)))
    saveas(gcf, strcat('C:\Users\derek\flow\trials\ke\feature\global-ke\s=0.05\', 'r=', string(fr), '.jpg'))
    close
end