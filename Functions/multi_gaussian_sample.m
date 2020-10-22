function sample = multi_gaussian_sample(centers, deviations, weights)
    len = length(weights);
    ind = randsample(1:len, 1, true, weights);
    sample = normrnd(centers(ind), deviations(ind));
end