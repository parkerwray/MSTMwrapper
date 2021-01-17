function fail = modified_kstest(radii, ropt, wopt, alpha)

[ropt, sort_idx] = sort(ropt);
wopt = wopt(sort_idx);

single_particle_wopt = zeros(length(wopt), 1);
single_particle_wopt(:) = wopt ./ (ropt .^ 2);
single_particle_wopt = single_particle_wopt ./ sum(single_particle_wopt);
%disp(single_particle_wopt)

cdf = zeros(length(ropt), 2);
cdf(:, 1) = ropt;
for idx = 1:length(ropt)
    cdf(idx, 2) = sum(single_particle_wopt(1:idx));
end
%figure,
%hist(radii, 1000);
%figure,
%bar(ropt, single_particle_wopt);




ref_radii = zeros(100000, 1);
ref_radii(1:(floor(single_particle_wopt(1) * length(ref_radii)))) = ropt(1);
%disp(floor(single_particle_wopt(1) * length(ref_radii)));
for i = 1:(length(single_particle_wopt) - 1)
    ref_radii((floor(cdf(i, 2) * length(ref_radii)) + 1):(floor(cdf(i+1, 2) * length(ref_radii)))) = ropt(i+1);
    %disp(length((floor(cdf(i, 2) * length(ref_radii)) + 1):(floor(cdf(i+1, 2) * length(ref_radii)))));
    %disp((floor(cdf(i, 2) * length(ref_radii)) + 1));
    %disp((floor(cdf(i, 2) * length(ref_radii)) + 1));
end
ref_radii((floor(cdf(end, 2) * length(ref_radii)) + 1):end) = ropt(end);
%{
ref_dist = @(~) randsample(ropt, 1, true, single_particle_wopt);
ref_radii = zeros(length(radii) * 100, 1);
for idx = 1:length(ref_radii)
    ref_radii(idx) = ref_dist(0);
end
%}

emp_cdf_r = min(ropt):0.1:(max(ropt)+1);
for i = 1:length(emp_cdf_r)
    emp_cdf(i) = sum(radii < emp_cdf_r(i)) / length(radii);
    ref_emp_cdf(i) = sum(ref_radii < emp_cdf_r(i)) / length(ref_radii);
end

figure,
hold on
plot(emp_cdf_r, emp_cdf);
plot(emp_cdf_r, ref_emp_cdf);
stairs(cdf(:,1), cdf(:,2));
legend("empirical", "empirical reference", "real");
hold off
%disp(cdf)
%cdf(:,2) = cat(1, [0], cdf(1:end-1, 2));

[fail, p] = kstest2(radii, ref_radii, 'Alpha', alpha);
fprintf("2 sample p=%d\n", p);
[fail, p, ksstat] = kstest(radii, 'CDF', cdf, 'Alpha', alpha);

fprintf("1 sample p=%d\n", p);
fprintf("1 sample test statistic = %d", ksstat);
fprintf("biggest cdf difference = %d", max(cdf(2:end,2) - cdf(1:(end-1),end)));
%disp(cdf)
end