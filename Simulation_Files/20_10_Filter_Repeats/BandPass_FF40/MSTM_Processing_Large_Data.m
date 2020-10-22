%% This code is used to process the results of 2D particle film simulations of random particle size distributions and random particle placement. The code extracts the modes of the center particle as a funciton of center particle size, cluster fill fraction, and wavelength. Statistics are generated beased on different unique distributions and polarization. The resulting statistics of the mode profile is then compared to the isolated particle case. 
% % 
% %% SPECIFY RELEVANT DIRECTORIES 

SAVE_FLAG = 1;

% Paths to import programs.
core = '/home/elipaul';
addpath(genpath(strcat(core,'/hypnos/Codes/MSTMwrapper'))); % Keep these the same
addpath(genpath(strcat(core,'/hypnos/Codes/Matlab_Functions')));  % Keep these the same
addpath(genpath(strcat(core,'/hypnos/Codes/randomparticles')));  % Keep these the same

% Specify the directory where you will place input and output files for
% communicating with the MSTM fortran program. Save files will also go to
% this location. 
parentdir = uigetdir(core,' Specify Simulation Directory'); % for Linux
%parentdir = uigetdir('$ME','Specify Save Directory'); % for NERSC

% Change directory to location where input/output files are saved
oldFolder = cd(parentdir);
%% MAKE SURE ALL FILES ARE VALID
result = zeros(3, Nsimulations);
FLAG = zeros(1, Nsimulations);

parfor (counter = 1:Nsimulations, 20)
    fname{counter} = strcat('mstm_', sprintf( '%03d',counter));
    [result(:,counter), FLAG(counter)] = export_mstm_intermediate(fname{counter});
end
%%
if sum(FLAG) ~= 0
    disp(fname(FLAG == 1));
end
figure,
plot(result(1,:));
figure,
plot(result(2,:));
figure,
plot(result(3,:));
%% PROCESS MODE COEFFICIENTS, FAR FIELDS, AND EFFICIENCIES 

clearvars A B Ax Bx A0 Ax0 B0 Bx0 
clearvars Ipar Iper Ipar0 Iper0
clearvars qe qa qsi qsd
clearvars size_param sphere_coffs_par sphere_coeffs_per
clearvars max_order theta

max_order = 10;
theta = [0, pi];

% Preallocate for speed in parfor
A = zeros(L1, L2, L3, L4, L5, max_order, max_order + 1);
B = zeros(L1, L2, L3, L4, L5, max_order, max_order + 1);
Ax = zeros(L1, L2, L3, L4, L5, max_order, max_order + 1);
Bx = zeros(L1, L2, L3, L4, L5, max_order, max_order + 1);
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

% Run parfor for speed of large dataset using linear vector notation 
parfor (counter = 1:L1*L2*L3*L4, 20)
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
            disp(fname);
        end
    end
        %continue;
        sphere_coeffs_par = squeeze(sphere_coeffs_data(1));
        sphere_coeffs_per = squeeze(sphere_coeffs_data(2));

        [A0(counter,:, :, :), B0(counter,:, :, :),...
            Ax0(counter, :, :, :), Bx0(counter, :, :, :), size_param] = ...
                    process_convert_modes_v2(max_order, sphere_coeffs_par.', sphere_coeffs_per.');
          
        [Ipar0(counter,:), Iper0(counter,:)] = make_S_v2(A0(counter,:, :, :), B0(counter,:, :, :),...
            Ax0(counter,:, :, :), Bx0(counter,:, :, :), theta);  
        
        [qe(counter), qa(counter), qsi(counter), qsd(counter)] = get_efficiencies(sphere_coeffs_par);
  
end
%% Reshape modes, far fields, and efficiencies

for counter = 1:L1*L2*L3*L4
   [idx1, idx2, idx3, idx4] = ind2sub([L1, L2, L3, L4], counter);
   A(idx1, idx2, idx3, idx4, :, :, :) = A0(counter,:,:,:);
   B(idx1, idx2, idx3, idx4, :, :, :) = B0(counter,:,:,:);
   Ax(idx1, idx2, idx3, idx4, :, :, :)= Ax0(counter,:,:,:);
   Bx(idx1, idx2, idx3, idx4, :, :, :)= Bx0(counter,:,:,:);
   Ipar(idx1, idx2, idx3, idx4, :) = Ipar0(counter, :);
   Iper(idx1, idx2, idx3, idx4, :) = Iper0(counter, :);
   qe(idx1, idx2, idx3, idx4) = qe(counter);
   qa(idx1, idx2, idx3, idx4) = qa(counter);
   qsi(idx1, idx2, idx3, idx4)= qsi(counter);
   qsd(idx1, idx2, idx3, idx4)= qsd(counter);
   
end
clearvars A0 B0 Ax0 Bx0 Iper0 Ipar0 sphere_coeffs_data sphere_coeffs_par sphere_coeffs_per
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