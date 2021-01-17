%%
%Note: Used NP_20_6
ropt = [ ...
  896.2125, ...
  820.9744, ...
  896.2126, ...
  884.4626, ...
  529.0447, ...
  592.4505, ...
  548.0243, ...
  592.4506, ...
  592.4506, ...
  801.8451, ...
  873.2437, ...
  592.4506, ...
  592.4505, ...
  699.9161, ...
  753.3129, ...
  652.9390];

wopt = [ ...
    0.1224, ...
    0.0506, ...
    0.1223, ...
    0.1259, ...
    0.0113, ...
    0.0476, ...
    0.0782, ...
    0.0475, ...
    0.0474, ...
    0.0441, ...
    0.0550, ...
    0.0474, ...
    0.0474, ...
    0.0140, ...
    0.0327, ...
    0.1063];




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
fill_fraction = [0.456285];
wavelengths = 4000:100:6000;
%wavelengths = linspace(500,900,81);
Ndistributions = 25; % Number of unique distributions to simulate. >=50 is found to be best but requires longer time. 
polarizations = [0,90];
r_mean = min(ropt);
index_info = {"file", '/home/elipaul/hypnos/Codes/Refractive Index Info/GaSb_Adachi.csv'};
cores = 20;
seed = 1205201441;

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

