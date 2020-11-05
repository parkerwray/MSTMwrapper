function [spheres, ff] = randomly_placed_normally_distributed_kerker_spheres(...
    type, dimension, ff_desired, center_radius, r_mean, r_sigma, bounds, giggles, loud, distr)

%{
    This function randomly places spheres in a 2D or 3D geomeery where the
    center particle is deterministically placed at the origin and has a
    deterministically (user set) radius.

    The dirstibution of particle radius is set to a normal distribution.
    The distribution of particle placement is set by a uniform
    distribution.

    This code has been tested for bugs and has multiple self consistency
    checks implemented at runtime. If an error is caught, the code goes
    into debug mode at that locaiton. 

    I think it is better to have isolated/modularized functions that
    generate spheres for a particular goal. These functions can be
    independently tested then ported to be used in the electromagnetic
    code. 

%} 

    
    
    
    
    
% test_PBCs; % Pass
% test_fcc; % Pass
% test_flips; % Pass
% test_make_random_v2; % Pass
if nargin < 10
   distr = @(~) random('normal', r_mean, r_sigma); 
end
scale = bounds; % Only for case of spherical boundary conditions
r = r_mean;
%sigma = r_sigma;
%distr = @(~) random('normal', r, sigma);
margin = 0.01;

tic;
if strcmp(type, "sphere") == 1
    [radii, ff, Nspheres] = get_radii_and_ff_in_sphere(scale, r, ...
        center_radius, ff_desired, distr, margin, dimension, loud);
    cords = full_randomize_in_sphere(radii, scale*r, giggles, dimension, loud);
elseif strcmp(type, "film") == 1
    if dimension == 3
        [cords, bounds, a] = make_fcc_3D(r, bounds, loud);
    else
        [cords, bounds, a] = make_fcc_2D(r, bounds, loud);
    end
    [radii, ff, Nspheres] = get_radii_and_ff(bounds, a,...
        center_radius, ff_desired, distr, margin, dimension, loud);
    [radii, cords] = full_randomize_v2(cords, radii, bounds.*a, ...
        giggles, dimension, ff_desired, margin, loud);
else
    disp("Invalid geometry requested."); 
    keyboard;
end

if loud
figure, 
plot_radii(radii); %pass
end

% WHAT IS THIS FUNCTION?? THE FLAG GOES TO 1?? - Parker
%has_intersections = check_intersection(cords, radii); %pass

% if has_intersections == 1
%     disp("Failure! Particles have intersections!")
%     keyboard;
% end
%     
    
%%
if loud
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
end
spheres = [radii, cords];



end