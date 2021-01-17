%%
clear
clc
%%
x = 1:20;
% Uniform Distribution:
% w = ones(10,1) * 0.1;
% Arbitrary:
w = rand(length(x),1) + 0.5;
w = w / sum(w);
dist = @(~) randsample(x, 1, true, w);



cdf = zeros(length(x), 2);
cdf(:,1) = x;
%cdf(:,2) = 0.1:0.1:1;
for i = 1:length(w)
    cdf(i,2) = sum(w(1:i));
end
%test = @(sample, cdf) kstest_stat(sample, cdf);
% chi squared test
test = @chi2gofstat_discrete;
%test = @(sample, cdf) chi2gofstat_discrete(sample, cdf);
sample_size = 200;
num_samples = 10000;


%%
statistics = make_test_statistic_distribution(dist, cdf, test, ...
    sample_size, num_samples);


%%
chi_x = 0:0.2:max(statistics);
figure,
hold on
histogram(statistics, 'Normalization', 'pdf')
plot(chi_x, chi2pdf(chi_x, length(x)-1));
%histogram(sqrt(sample_size) * statistics, 300, 'Normalization', 'probability')

%%
custom_tester = @(sample, alpha) (sum(statistics > test(sample, cdf)) / length(statistics) < alpha);
%passer = @(sample, alpha) disp(sum(statistics >= test(sample, cdf)) / length(statistics));% < (1 - alpha));

%% Test that it matches significance level in correct distributions
results = zeros(1000, 1);
alpha = 0.05;
for rep = 1:length(results)
    sample = make_sample(dist, sample_size);
    %tests(rep) = test(sample, cdf);
    %prop(rep) = sum(statistics > tests(rep)) / length(statistics);
    results(rep) = custom_tester(sample, alpha);
    %passer(sample, alpha);
end


disp(sum(results) / length(results));

%% Test against distributions skewed in a particular direction
bias_min = 0.9;
bias_max = 1.1;
low_bias = linspace(bias_min, bias_max, length(x));
low_bias_w = low_bias' .* w;
low_bias_dist = @(~) randsample(x, 1, true, low_bias_w);
high_bias = linspace(bias_max, bias_min, length(x));
high_bias_w = high_bias' .* w;
high_bias_dist = @(~) randsample(x, 1, true, high_bias_w);

resorted_w = w(randperm(length(w)));
resorted_dist = @(~) randsample(x, 1, true, resorted_w);

squared_w = w .^ 2;
squared_dist = @(~) randsample(x, 1, true, squared_w);

num_reps = 1000;
low_bias_res  = zeros(num_reps, 1);
high_bias_res = zeros(num_reps, 1);
resorted_res = zeros(num_reps, 1);
squared_res = zeros(num_reps, 1);
alpha = 0.05;

for rep = 1:num_reps
    low_sample  = make_sample(low_bias_dist , sample_size);
    high_sample = make_sample(high_bias_dist, sample_size);
    resorted_sample = make_sample(resorted_dist, sample_size);
    squared_sample = make_sample(squared_dist, sample_size);
    %tests(rep) = test(sample, cdf);
    %prop(rep) = sum(statistics > tests(rep)) / length(statistics);
    low_bias_res(rep) = custom_tester(low_sample, alpha);
    high_bias_res(rep) = custom_tester(high_sample, alpha);
    resorted_res(rep) = custom_tester(resorted_sample, alpha);
    squared_res(rep) = custom_tester(squared_sample, alpha);
    %passer(sample, alpha);
    %high_bias_dist_results(rep) = test(high_sample, cdf);
end

fprintf("Low bias results: %2f failure rate\n", sum(low_bias_res) / num_reps);
fprintf("High bias results: %2f failure rate\n", sum(high_bias_res) / num_reps);
fprintf("Resorted results: %2f failure rate\n", sum(resorted_res) / num_reps);
fprintf("Squared results: %2f failure rate\n", sum(squared_res) / num_reps);
%%
%cdf_eval([2,5.5], cdf)
mega_reg_sample = make_sample(dist, sample_size * num_reps);
mega_low_sample = make_sample(low_bias_dist, sample_size * num_reps);
reg_p_val = sum(statistics > (test(mega_reg_sample, cdf))) / length(statistics);
low_p_val = sum(statistics > (test(mega_low_sample, cdf))) / length(statistics);
%disp(p_val)
fprintf("Mega Regular p = %d\n", reg_p_val);
fprintf("Mega Biased p = %d\n", low_p_val);
%%
function sample = make_sample(dist, sample_size)
    sample = zeros(sample_size, 1);
    for idx = 1:sample_size
        sample(idx) = dist(0);
    end
end

function stat = kstest_stat(sample, cdf)
    [~, ~, stat] = kstest(sample, 'CDF', cdf);
end

function cdf = cdf_eval(x, cdf_in)
    idx = zeros(length(x), 1);
    for i = 1:length(x)
        idx(i) = find(cdf_in(:,1) - x(i) >= 0, 1);
        if x(i) ~= cdf_in(idx(i), 1)
            idx(i) = idx(i) - 1;
        end
    end
    cdf = cdf_in(idx, 2);
    %disp(cdf)
end

function stat = chi2gofstat_discrete(sample, cdf)
    [~, ~, stat_struct] = chi2gof(sample, 'Ctrs', cdf(:,1), 'CDF', {@cdf_eval, cdf});
    stat = stat_struct.chi2stat;
end