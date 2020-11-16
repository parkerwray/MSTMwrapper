function [A, Ax, B, Bx, sphere_coeffs_par0, sphere_coeffs_per0, spheres, sphere_coeffs_full, file_stats] = ...
    MSTM_Filter_Simulation_Step(...
    type, dimension, bounds, ...
    r_mean, ropt, wopt, fill_fraction, wavelengths, Ndistributions, ...
    index_info, convergence, cores, seed, ...
    function_directory,MSTMwrapper_directory, randomparticles_directory, ...
    mstm_location_serial, sim_directory)

%{
    index_info format: {type, optional params}
    For constant, it would look like {"constant", n, k}
    while for files, it would be {"file", filename}
%}

rng(seed);

index_type = index_info{1};
if strcmp(index_type, "constant")
    n_const = index_info{2};
    k_const = index_info{3};
elseif strcmp(index_type, "file")
    matfile = index_info{2};
end

center_radiis = ropt;

% This method only works for filters with given optimal radii (ropt) and
% AREA weights of those particles (wopt). These weights need to be rescaled
% to individual particle weights, i.e. the probability that a randomly
% selected particle will have a specific radius.
single_particle_wopt = wopt ./ (ropt.^2);
single_particle_wopt = single_particle_wopt ./ sum(single_particle_wopt);
dist = @(~) @(~) randsample(ropt, 1, true, single_particle_wopt);


% For clarity, leaving singular in input despite legacy use of plural
fill_fractions = fill_fraction;
% Needs a dummy value for legacy reasons
r_sigma = 0;


%core = core_directory;

%{
This code generates randomly positioned spheres in a 2D plane with radii
given by a set distribution. Material properties are assigned to each
sphere, then electromagnetic simulaitons (MSTM) are run and the results are
stored in a structure. 

This code is designed to approximate the properties of an infinite 2D
random particle film by simulating multiple finite sized, but sufficiently
large, simulations. Multiple layers of parallel processing is used to
minimize computation time. 
 
The simulation sweeps over:
1)  Center particle radius.
2) Cluster fill fraction (2D film).
3) Unique particle distributions of (1) and (2) for statistics. 
4) Wavelength/material properties.
5) Plane wave polarization at normal incidence.
%}

%clear;
%clc;

%SPECIFY IF YOU WANT TO SAVE THE RESULTS OF THE SIMULATION
%SAVE_FLAG = 1;

%SPECIFY RELEVANT DIRECTORIES 
% Paths to import programs.
%addpath(genpath(strcat(core,'/hypnos/Codes/MSTMwrapper'))); 
addpath(genpath(MSTMwrapper_directory));
%addpath(genpath(strcat(core,'/hypnos/Codes/Matlab_Functions')));  
addpath(genpath(function_directory));
%addpath(genpath(strcat(core,'/hypnos/Codes/randomparticles')));  
addpath(genpath(randomparticles_directory));

% Location of the fortran compiled run program for MSTM code.
%mstm_location = strcat(core, '/hypnos/Codes/MSTMwrapper/mstm_parallel_ubuntu.out');
%mstm_location_serial = strcat(core, ...
%    '/hypnos/Codes/MSTMwrapper/mstm_serial_ubuntu.exe');

% Specify the directory where you will place input and output files for
% communicating with the MSTM fortran program. Save files will also go to
% this location. 
%parentdir = uigetdir('/home/elipaul/hypnos',' Specify Save Directory'); % for Linux
parentdir = sim_directory;
%parentdir = uigetdir('$ME','Specify Save Directory'); % for NERSC

% Change directory to location where input/output files are saved
oldFolder = cd(parentdir);

%SPECIFY GENERAL PARAMETERS FOR MSTM RUN
%Set flags for MSTM run
mstm_flags = struct('write_sphere_data', 1,...
    'store_translation_matrix', 0,...
    'normalize_scattering_matrix',0,...
    'fixed_or_random_orientation',0,...
    'calculate_scattering_coefficients',1,...
    'track_nearfield_iterations',1,...
    'calculate_near_field',0,...  
    'calculate_t_matrix', 0,...   
    'azimuth_average_scattering_matrix',0);

%Define parameters for solution convergence


%{
convergence = struct(...
    'mie_epsilon', 10^-5,...
    'translation_epsilon', 10^-5,...
    'solution_epsilon', 10^-5,...
    'max_number_iterations', 1000,...
    'plane_wave_epsilon',10^-5,...
    't_matrix_convergence_epsilon',10^-5,...    
    'sm_number_processors', 1000,...
    'iterations_per_correction', 30);
%}
% NOTE: near_field_translation_distance is hard coded

%Define parameters for the input beam.






input_beam = struct(...
    'incident_or_target_frame', 0,... 
    'gaussian_beam_width', [],...
    'beam_type', 0,... % 0 = plane wave, 1 = Gaussian beam
    'incident_azimuth_angle_deg', 0,... % Alpha
    'incident_polar_angle_deg', 0,... % Beta 90 = supression
    'polarization_angle_deg',0); % Gamma

%Define parameters for near field calculations
near_field = struct(...
    'plane_cord',2,...
    'plane_position',0,...
    'plane_vertices',[-3000,-3000,3000,6000],...
    'resolution',0.5,...
    'near_field_output_data',2);

%Define parameters for particle distribution



mstm_input_params = struct(...
    'Nspheres', [],... 
    'k',[],... % length_scale_factor 
    'real_ref_index_scale_factor', 1,...
    'imag_ref_index_scale_factor',1,...
    'real_chiral_factor', 0,...
    'imag_chiral_factor', 0,...
    'medium_real_ref_index',1,...
    'medium_imag_ref_index',0,...
    'medium_real_chiral_factor',0,...
    'medium_imag_chiral_factor',0);

%SPECIFY SWEEP PARAMETERS FOR MSTM TO RUN
%Goal:
%Our goal is to look at the mode profiles in particles as they interact with other particles based on random radii and position. We then compare that result to the case if that  particle were isolated in free space. Our goal is to study how the interparticle coupling if effecting the solution and how we can use the results to predict the behavior of infinate (very large) films of these particles. 
%Motivation:
%Modern nm-sized paterning/lithography is not scalable to large area [> um^2] as drift/small inaccuracies in position compound to create large changes. Furthermore fabrication always has tolerance issues with respect to accuracy in shape. Metamaterials can be very sensitive to these facts if they rely on specific shape and periodicity to generate the desired effect. If particles can be found that reduce these dependencies, they could possibly be implemented over ~mm^2- m^2 area. This has particular applications in solar energy research. Other photonic/metamaterial applications apply as well. 
%Parameters in this code:
%We seek to compare the modes in the center particle to its isolated counterpart. Since we are comparing to the isolated case, we need to deterministically set the radius of the center particle so we can compare apples-to-apples. We then run many different random distributions around that center particle to get the statistics on how random coupling is effecting the result. These statistics not only tell us how the particle will respond to random coupling, on average, but have a direct link to a proposed formula for predicting an infinate film response. Currently there is no complete formula for defining the optical response in an infinate random film using the MSTM method (fundamentally it is likely not possible to define). We sweep the following:
%Center particle radius: since our film is composed of random sizes, we want to know how each particle size is effected by the random coupling with respect to its isolated particle counterpart. 
%Fill fraction: Coupling properties will change as a funciton of packing density. 
%Unique distributions: We need to run statistics for each center particle radius. This is achieved by simulating unique distributions around the center particle.
%Wavelength: We are interested in the modes and material properties as a function of wavelength. 
%Polarization: Averaging over the two polarization states dramatically simplifies the mathematics and is valid when talking about random films. 
%Distribution parameters:  Looking at different types of particle distributions (e.g., normal or bimodal). Studying this parameter should likely be done in seperate code as each new distribution is likely an independent study case. THIS CODE DOES NOT SWEEP DISTRIBUTION PARAMETERS.

% Good practice to always ensure variables are cleared before assigning. 
%clearvars r_mean r_sigma type dimension bounds giggles
%clearvars center_radiis fill_fractions Ndistributions wavelengths polarizations
clearvars k L1 L2 L3 L4 L5 Nsimulations

giggles = 500; % Putting this before import so it can be overwritten as desired
polarizations = [0, 90]; % Same as above (but it's a bad idea to change)


k = 2*pi./wavelengths;
L1 = length(center_radiis);
L2 = length(fill_fractions);
L3 = Ndistributions;
L4 = length(wavelengths);
L5 = length(polarizations);

Nsimulations = L1*L2*L3*L4*L5;
%disp(['Number of simulations: ', num2str(Nsimulations)])

%LOOP OVER PARAMETERS YOU WANT TO SWEEP

%NOTE: We are parallelizing on the number of unique distributions. This is usefull when we have >20 distributions, so we can make a maximum use of the parallel process. If not, we would likely want to parallelize on the wavelengths. This code sweeps a lot of parameters. If we do not parallelize these sweeps the code would take WAY too long. 

%Make the sphere distributions with material properties

%Run files here that generate particle distributions that you want to simulate. It is best that this code is determined by a secondary file that can be called as a function. The output must have the standard format for a unique distribution of spheres:
%Output:
%spheres = Nspheres x  [r, x, y, z, n, k] 

% You need to repeat static parameters 
mstm_input_params = repmat(mstm_input_params, L1, L2, L3, L4);
input_beam = repmat(input_beam, L5);

clearvars spheres ff
for idx1 = 1:L1 % Center radiis
    clearvars center_radii
    center_radii = center_radiis(idx1);
    %disp(strcat("center radius is ", num2str(center_radii)));
    for idx2 = 1:L2 % Fill fractions
        clearvars ff_desired
        ff_desired = fill_fractions(idx2);
        for idx3 = 1:L3 % Number of distributions
        %{
          Generate a particle distribution for a given desired fill fraction and center particle radius. Nspheres x [r, x, y, z]  
        %}
            clearvars shere_cords % Cord size may change, ensure variable is completely overwritten.  
            loud = 0;
            real_dist = dist(center_radii);
            [sphere_cords, ff(idx2, idx3)] = ...
                 randomly_placed_normally_distributed_kerker_spheres(...
                 type,...
                 dimension,...
                 ff_desired,...
                 center_radii,...
                 r_mean,...
                 r_sigma,...
                 bounds,...
                 giggles,...
                 loud,...
                 real_dist);
                  
            for idx4 = 1:L4 % Wavelengths
            %{
            Add wavelength dependent material properties. Save everything 
            to a cell based on the loop variables. Inside each
            cell element is [r x y z n k] for each N sheres. I.e., Nx6 
            %}
                %clearvars matfile
                if (strcmp(index_type, "file") == 1)
                    spheres(idx1, idx2, idx3, idx4).distribution = give_material_properties(sphere_cords,...
                                            matfile,...
                                            wavelengths(idx4));  
                elseif (strcmp(index_type, "constant") == 1)
                    spheres(idx1, idx2, idx3, idx4).distribution = [sphere_cords, ...
                        n_const.*ones(size(sphere_cords,1),1), k_const.* ones(size(sphere_cords,1),1)];
                else
                    disp("Invalid refractive index type!");
                    keyboard;
                end
                % Input params deals with lda and Nspheres so it need to
                % index both.
                mstm_input_params(idx1, idx2, idx3, idx4).k = k(idx4);       
                mstm_input_params(idx1, idx2, idx3, idx4).Nspheres = ...
                    size(spheres(idx1, idx2, idx3, idx4).distribution,1);
                                           
            end
        end
    end
end
    
for idx5 = 1:L5
    % Polarization sweep is not linked to other parameters. Remove from
    % the other loop to prevent unnecessary copying. 
    input_beam(idx5).incident_azimuth_angle_deg = polarizations(idx5);

end


%Parallelize the saving of files for MSTM to read.
%Data is written to text files to be read by the MSTM program. The writing process is parallelized using linear index notation. This allows for the maximum ammount of generality in the number of parameters you want to sweep as well as the maximum use of the parallelization cores. The only thing that need to be changed by the user are the output of the ind2sub function and the relevant indexes that are relevant to the sweep you desire. There is no need to add atitional for loops. 
 %disp('Generating input files for MSTM')
 %tic
parfor (counter = 1:Nsimulations,cores)
   
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
 %disp(['Files generated! Time: ', num2str(toc/60),'min'])
 %{
% Save setting used to make input files, if desired. 
if SAVE_FLAG == 1
    %Put letter "A" so the files are at the top of the folder. 
    save('A_Simulation_Input_Workspace.mat');
    save('A_Simulation_Input_File.m');
end
%}
%RUN SIMULATION FILES
%Uncomment if running on NERSC
%nodes = '2'; time = '00:05:00'; jobs = idx; make_mstm_SLURM_KNL_array_file(parentdir, nodes, time, jobs)
%Uncomment if running on Linux
% Copy mstm program file to folder
%copyfile(mstm_location, parentdir);
copyfile(mstm_location_serial, parentdir);

%disp('Running simulations.')
%tic
% Run parallelized simulation at a time
parfor (idx = 1:length(fname),cores)
    % Run MSTM code (Specify number of cores here "-n #")
    %command{idx} = ['/usr/lib64/openmpi/bin/mpirun -n 2 ./mstm_parallel_ubuntu.out ',fname{idx},'.inp'];
    command{idx} = ['./mstm_serial_ubuntu.exe ',fname{idx},'.inp'];
    %[status,cmdout] = system(command{idx},'-echo'); 
    [status,cmdout] = system(command{idx});
end
%disp(['Simulations finished! Time: ', num2str(toc/60),'min'])

%% MAKE SURE ALL RESULT FILES ARE VALID
result = zeros(3, Nsimulations);
parfor (counter = 1:Nsimulations, cores)
    fname{counter} = strcat('mstm_', sprintf( '%03d',counter));
    [result(:,counter), FLAG(counter)] = export_mstm_intermediate(fname{counter});
end
file_stats = result;
%%
if sum(FLAG) ~= 0
    throw(MException("Simulation:FileError", "Simulation output files are broken or missing!"));
end

%% PROCESS MODE COEFFICIENTS, FAR FIELDS, AND EFFICIENCIES 

clearvars A B Ax Bx A0 Ax0 B0 Bx0 
clearvars Ipar Iper Ipar0 Iper0
clearvars qe qa qsi qsd
clearvars size_param sphere_coffs_par sphere_coeffs_per
clearvars max_order theta

max_order = 10;
%theta = [0, pi];

% Preallocate for speed in parfor
A = zeros(L1, L2, L3, L4, L5, max_order, max_order + 1);
B = zeros(L1, L2, L3, L4, L5, max_order, max_order + 1);
Ax = zeros(L1, L2, L3, L4, L5, max_order, max_order + 1);
Bx = zeros(L1, L2, L3, L4, L5, max_order, max_order + 1);
%{
Iper = zeros(L1, L2, L3, L4, L5);
Ipar = zeros(L1, L2, L3, L4, L5);
A0 = zeros(L1*L2*L3*L4, L5, max_order, max_order+1);
Ax0 = zeros(L1*L2*L3*L4, L5, max_order, max_order+1);
B0 = zeros(L1*L2*L3*L4, L5, max_order, max_order+1);
Bx0 = zeros(L1*L2*L3*L4, L5, max_order, max_order+1);
Iper0 = zeros(L1*L2*L3*L4, L5);
Ipar0 = zeros(L1*L2*L3*L4, L5);
qe = zeros(L1, L2, L3, L4);
qa = zeros(L1, L2, L3, L4);
qsi = zeros(L1, L2, L3, L4);
qsd = zeros(L1, L2, L3, L4);
%}
% Run parfor for speed of large dataset using linear vector notation 
parfor (counter = 1:L1*L2*L3*L4, cores)
    %if mod(counter, 1000) == 0
    %    disp(counter);
        
    %end
    [idx1, idx2, idx3, idx4] =ind2sub([L1, L2, L3, L4], counter);
    
    sphere_coeffs_data = cell(2, 1);
    % Convert linear index notation to vector index notation 
    for idx5 = 1:L5
        
        fname =...
                strcat('mstm_',... 
                sprintf( '%03d', sub2ind([L1, L2, L3, L4, L5], idx1, idx2, idx3, idx4, idx5)));

        filename = strcat(parentdir,'/',fname);
        % DO NOT SAVE sphere_coeffs_data for all sweeps because the data is
        % too large. Instead process it and only save the processed
        % variables. 
        try
        [sphere_coeffs_data{idx5}] =...
            export_mstm_scattering_coeffs_v3(strcat(filename,'_scat_coeffs.dat'));
        catch
            throw(MException("Simulation:FileError", ...
                sprintf("File %s could not be read!", fname)));
        end
        sphere_coeffs_full{counter, idx5} = sphere_coeffs_data{idx5};
    end
        %continue;
        sphere_coeffs_par0{counter} = squeeze(sphere_coeffs_data(1));
        sphere_coeffs_per0{counter} = squeeze(sphere_coeffs_data(2));

        [A0(counter,:, :, :), B0(counter,:, :, :),...
            Ax0(counter, :, :, :), Bx0(counter, :, :, :), size_param] = ...
                    process_convert_modes_v2(max_order, sphere_coeffs_par0{counter}.', sphere_coeffs_per0{counter}.');
        %{ 
        [Ipar0(counter,:), Iper0(counter,:)] = make_S_v2(A0(counter,:, :, :), B0(counter,:, :, :),...
            Ax0(counter,:, :, :), Bx0(counter,:, :, :), theta);  
        
        [qe(counter), qa(counter), qsi(counter), qsd(counter)] = get_efficiencies(sphere_coeffs_par);
                    %}
  
end

%% Reshape modes, far fields, and efficiencies
sphere_coeffs_par = cell(L1, L2, L3, L4, L5);
sphere_coeffs_per = cell(L1, L2, L3, L4, L5);

for counter = 1:L1*L2*L3*L4
   [idx1, idx2, idx3, idx4] = ind2sub([L1, L2, L3, L4], counter);
   A(idx1, idx2, idx3, idx4, :, :, :) = A0(counter,:,:,:);
   B(idx1, idx2, idx3, idx4, :, :, :) = B0(counter,:,:,:);
   Ax(idx1, idx2, idx3, idx4, :, :, :)= Ax0(counter,:,:,:);
   Bx(idx1, idx2, idx3, idx4, :, :, :)= Bx0(counter,:,:,:);
   sphere_coeffs_par{idx1, idx2, idx3, idx4, 1} = sphere_coeffs_par0{...
       sub2ind([L1, L2, L3, L4, L5], [idx1, idx2, idx3, idx4, 1])};
   sphere_coeffs_par{idx1, idx2, idx3, idx4, 2} = sphere_coeffs_par0{...
       sub2ind([L1, L2, L3, L4, L5], [idx1, idx2, idx3, idx4, 1])};
   sphere_coeffs_per{idx1, idx2, idx3, idx4, 1} = sphere_coeffs_per0{...
       sub2ind([L1, L2, L3, L4, L5], [idx1, idx2, idx3, idx4, 1])};
   sphere_coeffs_per{idx1, idx2, idx3, idx4, 2} = sphere_coeffs_per0{...
       sub2ind([L1, L2, L3, L4, L5], [idx1, idx2, idx3, idx4, 1])};
   %{
   Ipar(idx1, idx2, idx3, idx4, :) = Ipar0(counter, :);
   Iper(idx1, idx2, idx3, idx4, :) = Iper0(counter, :);
   qe(idx1, idx2, idx3, idx4) = qe(counter);
   qa(idx1, idx2, idx3, idx4) = qa(counter);
   qsi(idx1, idx2, idx3, idx4)= qsi(counter);
   qsd(idx1, idx2, idx3, idx4)= qsd(counter);
   %}
end
clearvars A0 B0 Ax0 Bx0 Iper0 Ipar0 sphere_coeffs_data sphere_coeffs_par sphere_coeffs_per
A = squeeze(cat(dim_Ndistributions, A(:,:,:,:,1,:,:),A(:,:,:,:,2,:,:)));
B = squeeze(cat(dim_Ndistributions, B(:,:,:,:,1,:,:),B(:,:,:,:,2,:,:)));
Ax = squeeze(cat(dim_Ndistributions, Ax(:,:,:,:,1,:,:),Ax(:,:,:,:,2,:,:)));
Bx = squeeze(cat(dim_Ndistributions, Bx(:,:,:,:,1,:,:),Bx(:,:,:,:,2,:,:)));
%Legacy processing code:
%{
%% GET STATISTICS

clearvars dim_Ndistributions dim_wavelengths

dim_Ndistributions = 3;
dim_wavelengths = 4;

%% GET MODE STATISTICS

clearvars Astats Bstats Axstats Bxstats

% We want to treat the polarization data as an added number of
% distributions. I.e., a polarization of 90 degrees is the same as a
% distribution rotated by 90 degrees. Prior to this step, it was necessary
% to differentiate between a unique distribution and a flip in polarization
% to properly calculate the polarization conversion coefficients. Since
% these are calculated in the step above, we no longer need to
% differentiate. 
A = squeeze(cat(dim_Ndistributions, A(:,:,:,:,1,:,:),A(:,:,:,:,2,:,:)));
B = squeeze(cat(dim_Ndistributions, B(:,:,:,:,1,:,:),B(:,:,:,:,2,:,:)));
Ax = squeeze(cat(dim_Ndistributions, Ax(:,:,:,:,1,:,:),Ax(:,:,:,:,2,:,:)));
Bx = squeeze(cat(dim_Ndistributions, Bx(:,:,:,:,1,:,:),Bx(:,:,:,:,2,:,:)));

% Get average and standard deviation of mode energy and phase. Note: you
% need to first calculate the energy and phase then run the statistics,
% since energy and phase are nonlinear operators. 
[Astats] = get_mode_statistics(A, dim_Ndistributions);
[Bstats] = get_mode_statistics(B, dim_Ndistributions);
[Axstats] = get_mode_statistics(Ax, dim_Ndistributions);
[Bxstats] = get_mode_statistics(Bx, dim_Ndistributions);

% Grab the range min/max of the average mode energy in the wavelength
% dimension. We can use this to determine if a mode is negligable. 
Astats.energy = get_range(Astats.mean.mag,dim_wavelengths);
Bstats.energy = get_range(Bstats.mean.mag,dim_wavelengths);
Axstats.energy = get_range(Axstats.mean.mag,dim_wavelengths);
Bxstats.energy = get_range(Bxstats.mean.mag,dim_wavelengths);



%% GET PARTICLE EFFICIENCY STATISTICS

clearvars qestats qastats qsistats qsdstats

% Get average and standard deviation of particle efficiency.
[qestats] = get_statistics(qe, dim_Ndistributions);
[qastats] = get_statistics(qa, dim_Ndistributions);
[qsistats] = get_statistics(qsi, dim_Ndistributions);
[qsdstats] = get_statistics(qsd, dim_Ndistributions);


%% SAVE PROCESSED DATA

if SAVE_FLAG
    savedir = uigetdir(core,' Specify Save Directory'); % for Linux
    old_location = cd(savedir);
    mkdir('saved_data');
    cd('saved_data');
    save('modes.mat','A','B','Ax','Bx', '-v7.3');
    save('Is.mat', 'Ipar', 'Iper', '-v7.3');
    save('qs.mat', 'qe','qa','qsi','qsd','-v7.3');
    %save('sphere_distributions.mat', '-v7.3');
    save('simulation_parameters.mat', 'fill_fractions', 'center_radiis', 'Ndistributions', 'max_order', 'wavelengths', '-v7.3');    % UPDATE THIS LINE WITH THE RELEVANT PARAMETERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    save('modes_stats.mat','Astats','Bstats','Axstats','Bxstats','-v7.3');
    save('qs_stats.mat','qestats','qastats','qsistats','qsdstats','-v7.3');
end
%}
