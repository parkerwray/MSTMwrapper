


ropt = normrnd(100, 30, 100, 1);
wopt = rand(100, 1);
wopt = wopt .* (ropt .^ 2);
wopt = wopt ./ sum(wopt);
%{
wopt = 0.05:0.05:1;
ropt = norminv(wopt) + 20;
wopt = wopt .* (ropt .^ 2);
wopt = wopt ./ sum(wopt);
%}
single_particle_wopt = wopt ./ (ropt .^ 2);
single_particle_wopt = single_particle_wopt ./ sum(single_particle_wopt);
%single_particle_wopt = single_particle_wopt .^ 2;
%single_particle_wopt = single_particle_wopt ./ sum(single_particle_wopt);
%disp(single_particle_wopt);
dist = @(~) randsample(ropt, 1, true, single_particle_wopt);
num_repeats = 1;
results = zeros(num_repeats, 1);

alpha = 0.05;





[sorted_ropt, sort_idx] = sort(ropt);
sorted_s_wopt = single_particle_wopt(sort_idx);
%{
%plot(sorted_ropt, sorted_s_wopt);
figure,
%hold on
bar(ropt, single_particle_wopt);
figure,
bar(sorted_ropt, sorted_s_wopt);
%hold off

%}

%%

for rep = 1:length(results)
radii = zeros(10000, 1);
for idx = 1:length(radii)
    radii(idx) = dist(0);
end



%radii = (radii) / 2;

%radii = (radii + 100) / 2;
results(rep) = modified_kstest(sort(radii), ropt, wopt, alpha);

end
%{
figure,
hold on
plot(sorted_ropt, sorted_s_wopt, 'LineWidth', 6);
%ylim([0, 0.1]);
yyaxis right
histogram(radii, 1000, 'Normalization', 'probability');
%ylim([0, 0.1]);
hold off
%}

disp(sum(results) / length(results))