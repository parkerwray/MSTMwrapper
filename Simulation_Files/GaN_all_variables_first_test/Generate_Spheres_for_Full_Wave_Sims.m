
%{
Make particle distributions for full-wave lumerical simulations. 
%}

% Paths to import programs.
core = '/home/parkerwray';
addpath(genpath(strcat(core, '/hypnos/Codes/MSTMwrapper'))); % Keep these the same
addpath(genpath(strcat(core,'/hypnos/Codes/Matlab Functions')));  % Keep these the same
addpath(genpath(strcat(core,'/hypnos/Codes/randomparticles')));  % Keep these the same

% Good practice to always ensure variables are cleared before assigning. 
clearvars r_mean r_sigma type dimension bounds giggles

% Static parameters (use singular name conventiona):
r_mean = 140;
r_sigma = 10;
type = 'film';
dimension = 2;
bounds = [7, 7, 1];
giggles = 5000;
% Sweep parameters (user plural name convention):
center_radii = 140; % linspace(r_mean-r_sigma, r_mean+r_sigma, 3); % spacing set to have r = 140nm, to compare with  previous work
% CHOSE CENTER_RADII SPACING CAREFULLY TO MAXIMIZE USEFULLNESS IN STUDYING THE RESULT! 
ff_desired = 0.4;
loud = 1;

k = 2*pi./wavelengths;
clearvars spheres ff

[sphere_cords, ff] = ...
     randomly_placed_normally_distributed_kerker_spheres_v2(...
     type,...
     dimension,...
     ff_desired,...
     center_radii,...
     r_mean,...
     r_sigma,...
     bounds,...
     giggles,...
     loud);


    















