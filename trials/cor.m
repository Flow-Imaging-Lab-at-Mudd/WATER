function rsq = cor(est, real)

est = reshape(est, size(real));
    
rsq = 1 - sum((est - real).^2) / sum((real - mean(real)).^2);