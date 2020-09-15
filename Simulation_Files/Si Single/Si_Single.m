% The purpose of this simulation is to look at isolated Si particles over a
% wide variety of radii.

% This is the number of cores used for all computations. Make sure that
% matlab settings allow for at least this many to be used.
cores = 20;
% This is the radius on which everything is scaled--even if the
% distribution is not normal, it is important to specify this because 
% this is what the bounds are scaled based on. In general, the mean
% should be at least as small as the actual mean if not smaller. If errors
% are thrown about having insufficient particles, reducing mean and scaling 
% up bounds equivalently should fix the problem.
r_mean = 45;
% This does not matter except in the case of a normal distribution,
% but it is still a good idea to leave a dummy value for legacy reasons
r_sigma = 1;
% This specifies the type of region that is being simulated. The options
% are 'film' and 'sphere'. 'film' will give a cubic or square region, while
% 'sphere' will give a circular or spherical region.
type = 'film';
% This specifies the dimension in which the spheres are placed. The only
% currently implemented options are 2 and 3--see description of type for
% specifics.
dimension = 2;
% These are the bounds of the simulation region, given in lattice spacings
% based on r_mean. This means that the length will be r_mean*2sqrt(2)
% in both + and - directions.
bounds = [4, 4, 1];
% These are the center radiis which one wants to characterize. One should
% choose the spacing carefully to maximize usefullness--enough resolution
% to conclude about all relevant radii, but not so much that time is wasted
% in studying essentially equivalent cases.
center_radiis = 30:5:100;
% These are the fill fractions you want to study. This is defined as the
% n-dimensional cross-section which is occupied by spheres, where n is the dimension
% previously specified. This means that for the 2d case it is area
% fill-fraction, while for the 3d case it is volume fill. Some high fill
% fractions (>0.5 for 3D, >0.6 for 2D) may not be supported and will result
% in an infinite loop.
fill_fractions = [0, 0];
% These are the wavelengths which should be studied. Like center_radiis, it
% is important to choose an appropriate resolution that does not waste
% computing power. 
wavelengths = 300:5:800;
% This is the number of unique distributions to simulate. In general it is
% best to do at least 50, but for broad sweeps it can be advantageous to do
% fewer in the interest of time.
Ndistributions = 2;
% This is the radial distribution. It has two parts: the first input is the
% current radius, and the the second input is a null input, where it is a
% function that returns a random sample from the appropriate distribution.
% There are many possible approaches to defining a distribution--this is a
% very simple example, but the template for creating a distribution based
% on an area-weighted set of radii is given below this one.
dist = @(r) @(~) r;
%single_particle_wopt = wopt ./ (ropt.^2);
%single_particle_wopt = single_particle_wopt ./ sum(single_particle_wopt);
%dist = @(~) @(~) randsample(ropt, 1, true, single_particle_wopt);
%r_mean = sum(ropt .* single_particle_wopt)
%r_mean = min(ropt);

% This specifies how the refractive index is calculated. If 'file' is
% chosen, then matfile will be used in the calculation for each wavelength.
% If 'constant' is chosen instead, then n_const and k_const will be used in
% all cases instead.
index_type = 'file';
matfile = '/home/elipaul/hypnos/Codes/Refractive Index Info/Si_Aspnes.csv';
n_const = 4;
k_const = 0;


% Results: The kerker pattern continues to emerge over a huge range,
% although the backward region gets weaker for higher radii, while the
% forward region gets stronger. If the wavelengths get too low, higher
% order modes completely dominate.