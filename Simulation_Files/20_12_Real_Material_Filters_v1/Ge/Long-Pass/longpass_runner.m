%%
% Note: Used NP_20_27
ropt = [
  749.2249, ...
  565.8669, ...
  736.8252, ...
  712.1183, ...
  530.1508, ...
  502.0276, ...
  530.1508, ...
  791.2003, ...
  773.4692, ...
  478.2372, ...
  669.5459, ...
  530.1508, ...
  798.5562, ...
  502.0276, ...
  783.2299, ...
  699.5273, ...
  478.2371, ...
  530.1508, ...
  468.0189, ...
  893.6108];

wopt = [
    0.0599, ...
    0.0379, ...
    0.0578, ...
    0.0395, ...
    0.0410, ...
    0.0573, ...
    0.0409, ...
    0.1054, ...
    0.0268, ...
    0.0403, ...
    0.0437, ...
    0.0409, ...
    0.1428, ...
    0.0573, ...
    0.0630, ...
    0.0414, ...
    0.0404, ...
    0.0410, ...
    0.0228, ...
    0.0000];

%%
wopt_cutoff = 10^-4;
min_particle_diff = 1; %in nm
[ropt, wopt] = aggregate_similar_particles(ropt, wopt, wopt_cutoff, min_particle_diff);

%%
addpath(genpath('/home/elipaul/hypnos/Codes/MSTMwrapper')); 
MSTMwrapper_directory = '/home/elipaul/hypnos/Codes/MSTMwrapper';
%addpath(genpath('/home/elipaul/hypnos/Codes/Matlab_Functions'));  
function_directory = '/home/elipaul/hypnos/Codes/Matlab_Functions';
%addpath(genpath('/home/elipaul/hypnos/Codes/randomparticles'));  
randomparticles_directory = '/home/elipaul/hypnos/Codes/randomparticles';

% Location of the fortran compiled run program for MSTM code.
%mstm_location = '/home/elipaul/hypnos/Codes/MSTMwrapper/mstm_parallel_ubuntu.out';
mstm_location_serial = '/home/elipaul/hypnos/Codes/MSTMwrapper/mstm_serial_ubuntu.exe';

%%
% Specify the directory where you will place input and output files for
% communicating with the MSTM fortran program. Save files will also go to
% this location. 
sim_directory = uigetdir('/home/elipaul/hypnos',' Specify Save Directory'); % for Linux

convergence = struct(...
    'mie_epsilon', 10^-9,...
    'translation_epsilon', 10^-9,...
    'solution_epsilon', 10^-9,...
    'max_number_iterations', 1000,...
    'plane_wave_epsilon',10^-9,...
    't_matrix_convergence_epsilon',10^-9,...    
    'sm_number_processors', 1000,...
    'iterations_per_correction', 30);

%r_sigma = 100;
type = 'film';
dimension = 2;
bounds = [7, 7, 1];
%giggles = 250;
% Sweep parameters (user plural name convention):
center_radiis = ropt; % linspace(r_mean-r_sigma, r_mean+r_sigma, 3); % spacing set to have r = 140nm, to compare with  previous work
% CHOSE CENTER_RADII SPACING CAREFULLY TO MAXIMIZE USEFULLNESS IN STUDYING THE RESULT! 
fill_fraction = [0.443537];
wavelengths = 4000:100:10000;
%wavelengths = linspace(500,900,81);
Ndistributions = 25; % Number of unique distributions to simulate. >=50 is found to be best but requires longer time. 
polarizations = [0,90];
r_mean = min(ropt);
index_info = {"file", '/home/elipaul/hypnos/Codes/Refractive Index Info/Ge_Amotchkina.csv'};
cores = 20;
seed = 1221200025;

%%
% Note: file_stats was accidentally not included in the first run, thus is
% missing from the saved data
[A, Ax, B, Bx, sphere_coeffs_par0, sphere_coeffs_per0, spheres, sphere_coeffs_full, file_stats] = ...
    MSTM_Filter_Simulation_Step(...
    type, dimension, bounds, ...
    r_mean, ropt, wopt, fill_fraction, wavelengths, Ndistributions, ...
    index_info, convergence, cores, seed, ...
    function_directory,MSTMwrapper_directory, randomparticles_directory, ...
    mstm_location_serial, sim_directory);

