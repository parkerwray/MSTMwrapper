function spheres = example_sphere_function(type, dimension, ff_desired, r_mean, r_sigma)

%{
    I think it is better to have isolated/modularized functions that
    generate spheres for a particular goal. These functions can be
    independently tested then ported to be used in the electromagnetic
    code. 

%} 

% test_PBCs; % Pass
% test_fcc; % Pass
% test_flips; % Pass
% test_make_random_v2; % Pass
% dimension = 2;
% type = "sphere";
scale = 10;
r = r_mean;
sigma = r_sigma;
distr = @(~) random('normal', r, sigma);
% ff = 0.1;

margin = 0.01;

bounds = [2,2,2]; 
giggles = 100;
tic;
if strcmp(type, "sphere") == 1
    [radii, ff, Nspheres] = get_radii_and_ff_in_sphere(scale, r, ff_desired, ...
        distr, margin, dimension);
    cords = full_randomize_in_sphere(radii, scale*r, giggles, dimension);
elseif strcmp(type, "film") == 1
    if dimension == 3
        [cords, bounds, a] = make_fcc_3D(r, bounds);
    else
        [cords, bounds, a] = make_fcc_2D(r, bounds);
    end
    [radii, ff, Nspheres] = get_radii_and_ff(bounds, a,...
        ff_desired, distr, margin, dimension);
    [radii, cords] = full_randomize(cords, radii, bounds.*a, ...
        giggles, dimension);
else
    disp("Invalid geometry requested."); 
end
toc;
figure, 
plot_radii(radii); %pass
has_intersections = check_intersection(cords, radii); %pass
%%
% clc;
% lower_bound = bounds(1,:);
% upper_bound = bounds(2,:);
% Itop = find_top_view_ff(cords, r, lower_bound, upper_bound);
% fftop = sum(sum(Itop,1),2)/(size(Itop,1).*size(Itop,2));
% figure, imagesc(Itop)
% [Ntouch, ffcalc, is_zero, norm_dist, overlap_idx] = ...
%     check_touch_ff_orig(cords, ff, r, bounds(1,:), bounds(2,:));


%%
if strcmp(type, "film") == 1
    make_spheres(cords, radii, bounds(1,:).*a, bounds(2,:).*a);
    disp('Simulation region:')
    disp(['X distance: ', num2str((bounds(2,1)-bounds(1,1)).*a)])
    disp(['Y distance: ', num2str((bounds(2,2)-bounds(1,2)).*a)])
    disp(['Z distance: ', num2str((bounds(2,3)-bounds(1,3)).*a)])
elseif strcmp(type, "sphere") == 1
    make_spheres_in_sphere(cords, radii, r*scale);
    disp(['Cluster diameter: ',num2str((2*scale*r))])
    %make_spheres(cords, radii, [-1000,-1000,-1000], [1000,1000,1000]);
end

spheres = [radii, cords];



end