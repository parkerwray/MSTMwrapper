function [mean_val, std_val] = bootstrap_data(num_samples, N_tot, N_sample, ...
    A, Ax, B, Bx, Ipar, Iper, qa, qe, qsi, qsd, radii, wopt, ff, func)



for i = 1:num_samples
    data(i,:) = single_bootstrap_sample(N_tot, N_sample, A, Ax, B, Bx, ...
        Ipar, Iper, qa, qe, qsi, qsd, radii, wopt, ff, func);
end

mean_val = mean(data);
std_val = std(data);

end