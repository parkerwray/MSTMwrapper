function plot_deviations(data, dimension)

figure,
hold on
iterations = 1:size(data, dimension);
deviations = zeros(size(data, dimension), 1);
dat = squeeze(data);

for i = 1:length(iterations)
    deviations(i) = std(dat(1:i));
    means(i) = mean(dat(1:i));
end

yyaxis left
plot(iterations, deviations);
plot(iterations, std(dat) * ones(length(iterations), 1));
xlabel("Iterations");
ylabel("Standard Deviation");
yyaxis right
plot(iterations, means);
plot(iterations, mean(dat) * ones(length(iterations), 1));
ylabel("Mean")

end