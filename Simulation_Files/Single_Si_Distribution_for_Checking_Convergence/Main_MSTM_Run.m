%% This code generates randomly positioned spheres in a 2D plane with radii
%% given by a set distribution. Material properties are assigned to each
%% sphere, then electromagnetic simulaitons (MSTM) are run and the results are
%% stored in a structure. 
%% 
%% This code is designed to approximate the properties of an infinate 2D
%% random particle film by simulating multiple finite sized, but sufficiently
%% large, simulations. Multiple layers of parallel processing is used to
%% minimize computation time. 
%%  
%% The simulation sweeps over:
%% 1)  Center particle radius.
%% 2) Cluster fill fraction (2D film).
%% 3) Unique particle distributions of (1) and (2) for statistics. 
%% 4) Wavelength/material properties.
%% 5) Plane wave polarization at normal incidence.
%% 
%% 

clear;
clc;
%% SPECIFY IF YOU WANT TO SAVE THE RESULTS OF THE SIMULATION

SAVE_FLAG = 0;
%% SPECIFY RELEVANT DIRECTORIES 

% Paths to import programs.
core = '/home/parkerwray';
addpath(genpath(strcat(core, '/hypnos/Codes/MSTMwrapper'))); % Keep these the same
addpath(genpath(strcat(core,'/hypnos/Codes/Matlab Functions')));  % Keep these the same
addpath(genpath(strcat(core,'/hypnos/Codes/randomparticles')));  % Keep these the same


% Location of the fortran compiled run program for MSTM code.
mstm_location = strcat(core,'/hypnos/Codes/MSTMwrapper/mstm_parallel_ubuntu.out');
mstm_location_serial = strcat(core,'/hypnos/Codes/MSTMwrapper/mstm_serial_ubuntu.exe');

% Specify the directory where you will place input and output files for
% communicating with the MSTM fortran program. Save files will also go to
% this location. 
parentdir = uigetdir(core,' Specify Save Directory'); % for Linux
%parentdir = uigetdir('$ME','Specify Save Directory'); % for NERSC

% Change directory to location where input/output files are saved
oldFolder = cd(parentdir);
%% SPECIFY GENERAL PARAMETERS FOR MSTM RUN
%% Set flags for MSTM run

mstm_flags = struct('write_sphere_data', 1,...
    'store_translation_matrix', 0,...
    'normalize_scattering_matrix',0,...
    'fixed_or_random_orientation',0,...
    'calculate_scattering_coefficients',1,...
    'track_nearfield_iterations',1,...
    'calculate_near_field',0,...  
    'calculate_t_matrix', 0,...   
    'azimuth_average_scattering_matrix',0);
%% Define parameters for solution convergence

convergence = struct(...
    'mie_epsilon', 10^-9,...
    'translation_epsilon', 10^-9,...
    'solution_epsilon', 10^-9,...
    'max_number_iterations', 10000,...
    'plane_wave_epsilon',10^-9,...
    't_matrix_convergence_epsilon',10^-5,...    
    'sm_number_processors', 10000,...
    'iterations_per_correction', 30);
% NOTE: near_field_translation_distance is hard coded
%% Define parameters for the input beam.

% input_beam = struct(...
%     'incident_or_target_frame', 0,... 
%     'gaussian_beam_width', [],...
%     'beam_type', 0,... % 0 = plane wave, 1 = Gaussian beam
%     'incident_azimuth_angle_deg', 0,... % Alpha
%     'incident_polar_angle_deg', 0,... % Beta 90 = supression
%     'polarization_angle_deg',0); % Gamma
%% Define parameters for near field calculations

near_field = struct(...
    'plane_cord',2,...
    'plane_position',0,...
    'plane_vertices',[-3000,-3000,3000,6000],...
    'resolution',0.5,...
    'near_field_output_data',2);
%% Define parameters for particle distribution

% mstm_input_params = struct(...
%     'Nspheres', [],... 
%     'k',[],... % length_scale_factor 
%     'real_ref_index_scale_factor', 1,...
%     'imag_ref_index_scale_factor',1,...
%     'real_chiral_factor', 0,...
%     'imag_chiral_factor', 0,...
%     'medium_real_ref_index',1,...
%     'medium_imag_ref_index',0,...
%     'medium_real_chiral_factor',0,...
%     'medium_imag_chiral_factor',0);
%% SPECIFY SWEEP PARAMETERS FOR MSTM TO RUN

k = 2*pi./wavelengths;
L1 = 1; 
L2 = 1;
L3 = 1;
L4 = length(wavelengths);
L5 = length(polarizations);

Nsimulations = L1*L2*L3*L4*L5;
disp(['Number of simulations: ', num2str(Nsimulations)])

%% Make the sphere distributions with material properties

% 
% % You need to repeat static parameters 
% mstm_input_params = repmat(mstm_input_params, L1, L2, L3, L4);
% input_beam = repmat(input_beam, L5);
% 
% clearvars spheres ff
% for idx1 = 1:L1 % Center radiis
%     for idx2 = 1:L2 % Fill fractions
%         for idx3 = 1:L3 % Number of distributions
%             for idx4 = 1:L4 % Wavelengths
%                 mstm_input_params(idx1, idx2, idx3, idx4).k = k(idx4);       
%                 mstm_input_params(idx1, idx2, idx3, idx4).Nspheres = size(spheres(idx1, idx2, idx3, idx4).distribution,1);
%                                            
%             end
%         end
%     end
% end
%     
% for idx5 = 1:L5
%     % Polarization sweep is not linked to other parameters. Remove from
%     % the other loop to prevent unnecessary copying. 
%     input_beam(idx5).incident_azimuth_angle_deg = polarizations(idx5);
% 
% end

%% Parallelize the saving of files for MSTM to read.
% Data is written to text files to be read by the MSTM program. The writing 
% process is parallelized using linear index notation. This allows for the maximum 
% ammount of generality in the number of parameters you want to sweep as well 
% as the maximum use of the parallelization cores. The only thing that need to 
% be changed by the user are the output of the ind2sub function and the relevant 
% indexes that are relevant to the sweep you desire. There is no need to add atitional 
% for loops. 

 disp('Generating input files for MSTM')
 tic
parfor (counter = 1:Nsimulations,30)
   
    % Give a simulation file name based on linear index notation as a unique
    % identifier
    fname{counter} =...
            strcat('mstm_',... 
            sprintf( '%03d', counter ));
            
    % Convert linear index notation to vector index notation to send
    % relevant parameters to the input files for MSTM
    [idx1, idx2, idx3, idx4, idx5] =ind2sub([L1, L2, L3, L4, L5], counter);
        
    % This function generates MSTM inputs            
    make_mstm_job(parentdir,...
                        fname{counter},...
                        spheres(idx1, idx2, idx3, idx4).distribution,...
                        mstm_flags,...
                        convergence,...
                        input_beam(idx5),...
                        near_field,...
                        mstm_input_params(idx1, idx2, idx3, idx4)) 


end
 disp(['Files generated! Time: ', num2str(toc/60),'min'])
 
% Save setting used to make input files, if desired. 
if SAVE_FLAG == 1
    %Put letter "A" so the files are at the top of the folder. 
    save('A_Simulation_Input_Workspace.mat');
    save('A_Simulation_Input_File.mlx');
end
%% RUN SIMULATION FILES
%% Uncomment if running on NERSC

%nodes = '2'; time = '00:05:00'; jobs = idx; make_mstm_SLURM_KNL_array_file(parentdir, nodes, time, jobs)
%% Uncomment if running on Linux

% Copy mstm program file to folder
copyfile(mstm_location, parentdir);
copyfile(mstm_location_serial, parentdir);

disp('Running simulations.')
tic
% Run parallelized simulation at a time
parfor (idx = 1:length(fname),30)
    % Run MSTM code (Specify number of cores here "-n #")
    %command{idx} = ['/usr/lib64/openmpi/bin/mpirun -n 2 ./mstm_parallel_ubuntu.out ',fname{idx},'.inp'];
    command{idx} = ['./mstm_serial_ubuntu.exe ',fname{idx},'.inp'];
    %[status,cmdout] = system(command{idx},'-echo'); 
    if ~isfile([fname{idx}, '_output.dat'])
        
        [status,cmdout] = system(command{idx});
    end
end
disp(['Simulations finished! Time: ', num2str(toc/60),'min'])
%% EXTRACT MSTM DATA FOR POSTPROCESSING
% Data is read from the MSTM output text files and stored in relevant cells. 
% Thought the memory overhead is much larger using cells, it is necessary as the 
% output of each simulaiton is a structure with variables of different sizes inside. 
% The process of reading the MSTM output is based on parallelizing with linear 
% index notation. This allows for the maximum amount of generality on parameters 
% you want to sweep without needing to change the loop. The only necessary change 
% is the reshaping of the data.  

clearvars counter cluster_data0 sphere_data0 excitation_data0 sphere_coeff_data0
disp('Reading results.')
tic

parfor (counter = 1:length(fname),30)
        
    filename = strcat(parentdir,'/',fname{counter});
    try
%     [cluster_data0{counter},...
%         sphere_data0{counter},...
%         ~] =...
%         export_output_file_v2(strcat(filename,'_output.dat'));
    
    [sphere_coeffs_data0{counter}] =...
        export_mstm_scattering_coeffs_v3(strcat(filename,'_scat_coeffs.dat'));
    catch
        disp(counter)
    end

end
disp(['Results Read! Time: ', num2str(toc/60),'min'])


%% Reshape the data to match the relevant sweep parameters

clearvars cluster_data sphere_data sphere_coeffs
for counter = 1:length(fname)
    
    % Convert linear index notation to vector index notation to send
    % relevant parameters to the input files for MSTM
    [idx1, idx2, idx3, idx4, idx5] =ind2sub([L1, L2, L3, L4, L5], counter);
    %cluster_data{idx1, idx2, idx3, idx4, idx5} = cluster_data0{counter};
    %sphere_data{idx1, idx2, idx3, idx4, idx5} = sphere_data0{counter};
    sphere_coeffs_data{idx1, idx2, idx3, idx4, idx5} = sphere_coeffs_data0{counter};
    
end
clearvars cluster_data0 sphere_data0 sphere_coeffs_data0
%%
oldFolder = cd(parentdir);
if SAVE_FLAG == 1
    %save('A_Simulation_Output_Workspace.mat'); % VERY LARGE!
    save('A_simulation_coeffs.mat', "sphere_coeffs_data"); % To not have to open large file later.
end