


addpath(genpath('/home/elipaul/hypnos/Codes/MSTMwrapper')); 
addpath(genpath('/home/elipaul/hypnos/Codes/Matlab_Functions'));  
addpath(genpath('/home/elipaul/hypnos/Codes/randomparticles'));  

type = "film";
dimension = 2;
%resolution = 301;
Nsimulations = 200;
ff_goal = 0.4;
center_radius = 500;
bounds = [6, 6, 1];
r_mean = 500;
r_sigma = 50;
giggles = 100;
loud = 0;
dist = @(~) normrnd(r_mean, r_sigma);
%dist = @(~) multi_gaussian_sample([100, 250, 371, 500, 850, 1000], [10, 40, 80, 5, 200, 150], [0.1, 0.05, 0.1, 0.2, 0.3, 0.25]);
%dist = @(~) 1000 * rand(1,1);


half_width = bounds(1) * r_mean * 2 * sqrt(2);
resolution = floor(half_width * 2 / 125 + 1);

prob_dist = zeros(resolution, resolution);


total_spheres = 0;

for i = 1:Nsimulations
    if mod(i, 10) == 0
        disp(i)
    end
    [spheres, ~] = randomly_placed_normally_distributed_kerker_spheres(...
        type, dimension, ff_goal, center_radius, r_mean, r_sigma, bounds, ...
        giggles, loud, dist);
    for j = 1:length(spheres)
        % Convert from absolute coordinates to index coordinates
        x_ind = ceil(((spheres(j,2) / half_width) + 1) * resolution / 2);
        y_ind = ceil(((spheres(j,3) / half_width) + 1) * resolution / 2);
        
        % Ignore mirror particles
        if x_ind > 0 && y_ind > 0 && x_ind <= resolution && y_ind <= resolution && ~(spheres(j, 2) == 0 && spheres(j, 3) == 0)
            prob_dist(x_ind, y_ind) = prob_dist(x_ind, y_ind) + 1;
            
            total_spheres = total_spheres + 1;
        end
    end
        
end

prob_dist = prob_dist ./ total_spheres;


transformed = fft2(prob_dist);

%%
sample_arr = zeros(10000, 1);
for i = 1:10000
sample_arr(i) = dist(0);
end
figure,
histogram(sample_arr)
%%
center_ind = ceil(resolution / 2);
figure,
hold on
imagesc(prob_dist)
viscircles([center_ind, center_ind], (center_radius + min(sample_arr)) / 125)
%imagesc(log((total_spheres * prob_dist)+1))
colorbar
hold off
%figure,
%imagesc(log(abs(transformed)))
%colorbar



% Make dists and reshape


%%
