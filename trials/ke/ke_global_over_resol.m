% Spacings of velocity.
sps = 0.01: 0.05: 1;

for sp = sps
    KE_uniform
    saveas(gcf, strcat('s=', string(sp), '.jpg'))
    close
end