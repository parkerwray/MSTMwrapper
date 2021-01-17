function statistics = make_test_statistic_distribution(dist, cdf, test, ...
    sample_size, num_samples)

%{
    Input formats:
    dist: lambda function returning a random sample from the distribution
        when given an input of 0
    cdf: actual cdf of the distribution, in two-column form
    test: lambda function taking in (sample, cdf) and returning a test
        statistic
    sample_size: integer number of elements per sample
    num_samples: integer number of samples to repeat
%}

statistics = zeros(num_samples, 1);
for i = 1:num_samples
    statistics(i) = test(make_sample(dist, sample_size), cdf);
    %sample = make_sample(dist, sample_size);
    %[~, ~, statistics(i)] = test(sample, cdf);
end



function sample = make_sample(dist, sample_size)
    sample = zeros(sample_size, 1);
    for idx = 1:sample_size
        sample(idx) = dist(0);
    end
end

end