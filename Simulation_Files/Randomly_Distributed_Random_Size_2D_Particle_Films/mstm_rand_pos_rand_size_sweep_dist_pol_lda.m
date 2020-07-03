clear;
clc;

%{

This code generates randomly positioned spheres in a 2D plane with radii
given by a set distribution. Material properties are assigned to each
sphere, then electromagnetic simulaitons (MSTM) are run and the results are
stored in a structure. 

This code is designed to approximate the properties of an infinate 2D
random particle film by simulating multiple finite sized, but sufficiently
large, simulations. Multiple layers of parallel processing is used to
minimize computation time. 
 

The simulation sweeps over:
1) Unique particle distribution. 
2) Wavelength.
3) Polarization.

Inputs are:
1) Plane wave at normal incidence.


%}



%% GENERAL FLAGS
SAVE_FLAG = 0;

%% SPECIFY LOCATION TO SAVE INPUT FILES
%% Use if running on NERSC
% parentdir = uigetdir('$ME','Specify Save Directory'); %Use if running on NERSC cluster

%% Use if running on Windows
%parentdir = uigetdir('Z:','Specify Save Directory');
 
%% Use if running on Linux
addpath(genpath('/home/parkerwray/hypnos/Codes/MSTMwrapper'));  
addpath(genpath('/home/parkerwray/hypnos/Codes/Matlab Functions'));  
parentdir = uigetdir('/home/parkerwray/hypnos',' Specify Save Directory');

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
    'max_number_iterations', 100000,...
    'plane_wave_epsilon',10^-9,...
    't_matrix_convergence_epsilon',10^-9,...    
    'sm_number_processors', 1000,...
    'iterations_per_correction', 30);

%% Define parameters for the input beam.
input_beam = struct(...
    'incident_or_target_frame', 0,... 
    'gaussian_beam_width', [],...
    'beam_type', 0,... % 0 = plane wave, 1 = Gaussian beam
    'incident_azimuth_angle_deg', 0,... % Alpha
    'incident_polar_angle_deg', 0,... % Beta 90 = supression
    'polarization_angle_deg',0); % Gamma

%% Define parameters for near field calculations
near_field = struct(...
    'plane_cord',2,...
    'plane_position',0,...
    'plane_vertices',[-3000,-3000,3000,6000],...
    'resolution',0.5,...
    'near_field_output_data',2);

%% Define parameters for particle distribution
mstm_input_params = struct(...
    'Nspheres', [],... 
    'k',[],...
    'real_ref_index_scale_factor', 1,...
    'imag_ref_index_scale_factor',1,...
    'real_chiral_factor', 0,...
    'imag_chiral_factor', 0,...
    'medium_real_ref_index',1,...
    'medium_imag_ref_index',0,...
    'medium_real_chiral_factor',0,...
    'medium_imag_chiral_factor',0);


%% SPECIFY SWEEP PARAMETERS FOR MSTM TO RUN
%%
clearvars wavelengths Ndistributions L1 L2 L3

wavelengths = linspace(300,600,300); %2:0.05:10;  % len = 11
k = 2*pi./wavelengths;
Ndistributions = 2;
polarization = [0,90];

L1 = Ndistributions; % Number of unique distributions to simulate
L2 = length(wavelengths);
L3 = length(polarization);

disp(['Number of sims: ', num2str(L1*L2*L3)])


%% LOOP OVER PARAMETERS YOU WANT TO SWEEP
%%
%{
    NOTE: We are parallelizing on the number of unique distributions. This
    is usefull when we have >20 distributions, so we can make a maximum use
    of the parallel process. If not, we would likely want to parallelize on
    the wavelengths. 

%}



%% Construct file names based on a unique identifier
counter = 1;
for idx1 = 1:L1
    for idx2 = 1:L2
        for idx3 = 1:L3
           fname{counter} =...
                    strcat('mstm_',... %GIVE A FILE NAME!
                    sprintf( '%03d', counter ));
           
           counter_map(idx1, idx2, idx3) = counter;          
           counter = counter+1;     
            
        end
    end
end

%% Make the sphere distributions with material properties
clearvars spheres
for idx1 = 1:L1
    
   %{ 
        Run files here that generate particle distributions that you want to
        simulate. It is best that this code is determined by a secondary file that
        can be called as a function. The output must have the standard format:
        Input: 
            wavelengths
        Output:
            spheres = [r, x, y, z]
    %}
    
     type = 'film';
     dimension = 2;
     ff_desired = 0.2;
     r_mean = 44;
     r_sigma = 0.1*r_mean;
    
    sphere_cords = example_sphere_function(type,...
                          dimension,...
                          ff_desired,...
                          r_mean,...
                          r_sigma);
    
    for idx2 = 1:L2
        %{
            Add wavelength dependent material properties. Save everything 
            to a cell based on the loop variables. Inside each
            cell element is [r x y z n k] for each N sheres. I.e., Nx6 
        %}
        matfile = '/home/parkerwray/hypnos/Project - Dusty Plasma MURI/Material Data/Refractive Index Info/Si_Aspnes.csv';
        spheres{idx1, idx2} = give_material_properties(sphere_cords,...
                                           matfile,...
                                           wavelengths(idx2));                                 
    end                   
end


%% Parallelize the saving of files for MSTM to read. 
 disp('Generating input files for MSTM')
 tic
parfor (idx1 = 1:L1,20)
    for idx2= 1:L2
        
        for idx3 = 1:L3
            
            % Change relevant simulation settings for sweep
            dummy_input_beam = input_beam;
            dummy_mstm_input_params = mstm_input_params;
            
            dummy_input_beam.incident_azimuth_angle_deg = polarization(idx3);
            dummy_mstm_input_params.Nspheres = size(spheres{idx1, idx2},1);
            dummy_mstm_input_params.k = k(idx2);

            
            
            % This function generates MSTM inputs
            counter = counter_map(idx1, idx2,idx3);              
            make_mstm_job(parentdir,...
                                fname{counter},...
                                spheres{idx1, idx2},...
                                mstm_flags,...
                                convergence,...
                                dummy_input_beam,...
                                near_field,...
                                dummy_mstm_input_params) 

        end
    end
end
 disp(['Files generated! Time: ', num2str(toc/60),'min'])       

 
%% RUN SIMULATION FILES
%%

% Change directory to location where input files are saved
oldFolder = cd(parentdir);

% Save setting used to make input files, if desired. 
if SAVE_FLAG == 1
    save('Simulation_Input_Workspace.mat');
    save('Simulation_Input_File.m');
end

%% Use if running on NERSC
% nodes = '2';
% time = '00:05:00';
% jobs = idx;
% make_mstm_SLURM_KNL_array_file(parentdir, nodes, time, jobs)

%% Use if running on Linux

% Copy mstm program file to folder
mstm_location = '/home/parkerwray/hypnos/Codes/MSTMwrapper/mstm_parallel_ubuntu.out';
copyfile(mstm_location, parentdir);

disp('Running simulations.')
tic
% Run parallelized simulation at a time
parfor (idx = 1:length(fname),5)
    % Run MSTM code (Specify number of cores here "-n #")
    command{idx} = ['/usr/lib64/openmpi/bin/mpirun -n 3 ./mstm_parallel_ubuntu.out ',fname{idx},'.inp'];
    %[status,cmdout] = system(command{idx},'-echo'); 
    [status,cmdout] = system(command{idx});
end
disp(['Simulations finished! Time: ', num2str(toc/60),'min'])  

%% EXTRACT MSTM DATA FOR POSTPROCESSING
%%
% Parallelize data reading process
clearvars counter cluster_data sphere_data excitation_data
disp('Reading results.')
tic
parfor (idx1 = 1:L1,20)
    for idx2= 1:L2
        for idx3 = 1:L3
            counter = counter_map(idx1, idx2, idx3) ;   
            filename = strcat(parentdir,'/',fname{counter});

            [cluster_data{idx1, idx2, idx3},...
                sphere_data{idx1, idx2, idx3},...
                excitation_data{idx1, idx2, idx3}] =...
                export_output_file_v2(strcat(filename,'_output.dat'));
            
            [sphere_coeffs{idx1, idx2, idx3}] =...
                export_mstm_scattering_coeffs_v3(strcat(filename,'_scat_coeffs.dat'));


        end
    end
end
disp(['Results Read! Time: ', num2str(toc/60),'min'])  


%%
oldFolder = cd(parentdir);
if SAVE_FLAG == 1
    save('Simulation_Output_Workspace.mat');
end



% 
% 
% 
% %% Get absorption Efficiency
% clear qa
% for idx1 = 1:L1
%     for idx2= 1:L2
%        for idx3 = 1:L3
%         qa(idx1, idx2, idx3) =...
%             cluster_data{idx1, idx2, idx3}.qabs;
%        end
%     end
% end
% 
% 
% 
% %%
% 
% % mask = (qa(:,:,1) < 0.05);
% q_ratio = qa(:,:,1)./qa(:,:,2);
% % mask2 = (qa(:,:,2) < 0.25);
% q_ratio2 = qa(:,:,2)./qa(:,:,1);
% q_fc = (qa(:,:,1)-qa(:,:,2))./(qa(:,:,1)+qa(:,:,2));
% % q_ratio(mask) = NaN;
% % q_ratio2(mask2) = NaN;
% 
% 
% figure;        
% k = 2*pi./lda;
% rcore = xcore./k;
% us = 3;
% q_ratio = imgaussfilt(imresize(q_ratio,us),1.2);
% y = linspace(1.3,3,80*us);
% x = linspace(2,10,80*us); %2:0.05:10;  % len = 11
% imagesc(x, y, q_ratio)
% xlabel('x shell')
% ylabel('n shell')
% title(['x_{core} = ', num2str(xcore),newline,'m_{core} = ',num2str(ncore),'+',num2str(kcore),'i'])
% colormap parula
% h2 = colorbar;
% ylabel(h2, 'Forward-to-backward absorption/emission ratio')
% pbaspect([1 1 1])
% 
% 
% %%
% 
% figure;
% subplot(2,2,1)
% imagesc(xshell, nshell, q_ratio)
% xlabel('x shell')
% ylabel('n shell')
% title(['n core = ',num2str(ncore),'+',num2str(kcore),'i, n shell = ',num2str(nshell(IDX)),'+',num2str(kshell),'i'])
% colormap jet
% h2 = colorbar;
% ylabel(h2, '(0)/(180)')
% 
% subplot(2,2,2)
% imagesc(xshell, nshell, q_ratio2)
% xlabel('x shell')
% ylabel('n shell')
% title(['n core = ',num2str(ncore),'+',num2str(kcore),'i, n shell = ',num2str(nshell(IDX)),'+',num2str(kshell),'i'])
% colormap jet
% h2 = colorbar;
% ylabel(h2, '(180)/(0)')
% 
% subplot(2,2,3)
% imagesc(xshell, nshell, qa(:,:,1))
% xlabel('x shell')
% ylabel('n shell')
% title(['n core = ',num2str(ncore),'+',num2str(kcore),'i, n shell = ',num2str(nshell(IDX)),'+',num2str(kshell),'i'])
% colormap jet
% h2 = colorbar;
% ylabel(h2, 'Qa at theta = 0')
% 
% subplot(2,2,4)
% imagesc(xshell, nshell, qa(:,:,2))
% xlabel('x shell')
% ylabel('n shell')
% title(['n core = ',num2str(ncore),'+',num2str(kcore),'i, n shell = ',num2str(nshell(IDX)),'+',num2str(kshell),'i'])
% colormap jet
% h2 = colorbar;
% ylabel(h2, 'Qa at theta = 90')
% 
% % subplot(2,2,4)
% % imagesc(xshell, xcore(15:end), qa(15:end,:,3))
% % xlabel('x shell')
% % ylabel('x core')
% % title(['n core = 2+2i, n shell = 2+0i'])
% % h2 = colorbar;
% % ylabel(h2, 'Qa at theta = 180')
% %%
% % savefig(['Assym_Emission_Particle_Outside', num2str(COUNTER),'.fig'])
% % f = gcf;
% % 
% % saveas(f,['Assym_Emission_Particle_Outside', num2str(COUNTER),'.tif']);
% % %exportgraphics(f,['Assym_Emission', num2str(COUNTER),'.tif'],'Resolution',300)
% % COUNTER = COUNTER+1;
% %%
% 
% % [val, idx] = max(q_fc,[],2);
% % [val2, idx2 ] = max(val);
% % 
% % % [val, idx] = min(abs(22000-rshell));
% % figure,
% % hold on 
% % plot(xshell, squeeze(qa(1,:,1)))
% % plot(xshell, squeeze(qa(1,:,2)))
% % hold off
% % xlabel('xshell')
% % ylabel('Absorption efficiency')
% % title(['x core: ',...
% %     num2str(xcore(1))]);
% % legend('Qa at 0 deg', 'Qa at 180 deg');
% 
% 
% %%
% 
% % 
% figure,
% hold on 
% plot(wavelengths/1000, squeeze(qa(:,4,1)))
% plot(wavelengths/1000, squeeze(qa(:,4,2)))
% hold off
% xlabel('Wavelength (um)')
% ylabel('Absorption efficiency')
% title(['Core Material: SiO2 Core Radius: ',...
%     num2str(rcore/1000), 'um ', 'Shell Radius: ',...
%     num2str(rshell(4)/1000), 'um'])
% legend('Qa at 0 deg', 'Qa at 180 deg')





% %%
% oldFolder = cd(parentdir);
% if SAVE_FLAG == 1
%     save('Simulation_Output_Workspace.mat');
% end
% 



% %[status,cmdout] = system(command,'-echo');
% [status,cmdout] = system(command{idx});
% 
% 
% 
% 
% 
% % Run multiple parallelized simulation at a time
% commands = [];
% for idx = 1:length(fname)
%     % Run MSTM code (Specify number of cores here)
%     command = ['/usr/lib64/openmpi/bin/mpirun -n 10 ./mstm_parallel_ubuntu.out ',fname{idx},'.inp'];
%     if idx > 1
%     dummy = [commands, ' && ',command];
%         if mod(idx,
%     %[status,cmdout] = system(command,'-echo');
%     [status,cmdout] = system(command);
% end
% 
% 
% 