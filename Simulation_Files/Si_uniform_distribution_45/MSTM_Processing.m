%% This code is used to process the results of 2D particle film simulations of random particle size distributions and random particle placement. The code extracts the modes of the center particle as a funciton of center particle size, cluster fill fraction, and wavelength. Statistics are generated beased on different unique distributions and polarization. The resulting statistics of the mode profile is then compared to the isolated particle case. 
% % 
% %% SPECIFY RELEVANT DIRECTORIES 

core = '/home/elipaul';
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

%% PROCESS MODE COEFFICIENTS 


clearvars A B Ax Bx size_param sphere_coffs_par sphere_coeffs_per Ipar Iper
max_order = 10;
theta = [0, pi];
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
parfor (counter = 1:L1*L2*L3*L4, 10)
    if mod(counter, 1000) == 0
        disp(counter);
        
    end
    [idx1, idx2, idx3, idx4] =ind2sub([L1, L2, L3, L4], counter);
    
    sphere_coeffs_data = cell(2, 1);
    % Convert linear index notation to vector index notation 
    for idx5 = 1:L5
        
        fname =...
                strcat('mstm_',... 
                sprintf( '%03d', sub2ind([L1, L2, L3, L4, L5], idx1, idx2, idx3, idx4, idx5)));

        filename = strcat(parentdir,'/',fname);
        [sphere_coeffs_data{idx5}] =...
            export_mstm_scattering_coeffs_v3(strcat(filename,'_scat_coeffs.dat'));
    end

        sphere_coeffs_par = squeeze(sphere_coeffs_data(1));
        sphere_coeffs_per = squeeze(sphere_coeffs_data(2));
        
        [A0(counter,:, :, :), B0(counter,:, :, :),...
            Ax0(counter, :, :, :), Bx0(counter, :, :, :), size_param] = ...
                    process_convert_modes_v2(max_order, sphere_coeffs_par.', sphere_coeffs_per.');
          
        [Ipar0(counter,:), Iper0(counter,:)] = make_S_v2(A0(counter,:, :, :), B0(counter,:, :, :),...
            Ax0(counter,:, :, :), Bx0(counter,:, :, :), theta);  
        
        [qe(counter), qa(counter), qsi(counter), qsd(counter)] = get_efficiencies(sphere_coeffs_par);
  
end
for counter = 1:L1*L2*L3*L4
   [idx1, idx2, idx3, idx4] = sub2ind([L1, L2, L3, L4], counter);
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
clearvars A0 B0 Ax0 Bx0 Iper0 Ipar0

%% GET MODE STATISTICS


clearvars A_mag_stats A_phase_stats
clearvars B_mag_stats B_phase_stats
clearvars Ax_mag_stats Ax_phase_stats
clearvars Bx_mag_stats Bx_phase_stats

% We want to treat the polarization data as an added number of
% distributions. I.e., a polarization of 90 degrees is the same as a
% distribution rotated by 90 degrees. Prior to this step, it was necessary
% to differentiate between a unique distribution and a flip in polarization
% to properly calculate the polarization conversion coefficients. Since
% these are calculated in the step above, we no longer need to
% differentiate. 

dim_Ndistributions = 3;
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
dim_wavelengths = 4;
Astats.energy = get_range(Astats.mean.mag,dim_wavelengths);
Bstats.energy = get_range(Bstats.mean.mag,dim_wavelengths);
Axstats.energy = get_range(Axstats.mean.mag,dim_wavelengths);
Bxstats.energy = get_range(Bxstats.mean.mag,dim_wavelengths);

% Get average and standard deviation of particle efficiency.
[qestats] = get_statistics(qe, dim_Ndistributions);
[qatats] = get_statistics(qa, dim_Ndistributions);
[qsistats] = get_statistics(qsi, dim_Ndistributions);
[qsdstats] = get_statistics(qsd, dim_Ndistributions);


%% Save files
%Save As
%Save Is
%Save qs
%save spheres.distribution
%Save basic params
%Save stats

old_location = cd(parentdir);
mkdir('A_saved_data');
save('modes.mat','A','B','Ax','Bx', '-v7.3');
save('Is.mat', 'Ipar', 'Iper', '-v7.3');
save('qs.mat', 'qe','qa','qsi','qsd','-v7.3');
save('sphere_distributions','-v7.3');
save('simulation_parameters', 'ff','-v7.3');    % UPDATE THIS LINE WITH THE RELEVANT PARAMETERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
save('modes_stats.mat','Astats','Bstats','Axstats','Bxstats','-v7.3');
save('qs_stats.mat','qestats','qastats','qsistats','qsdstats','-v7.3');

















%% PLOT MAGNITUDES AND PHASES
% % 
% % 
% 

x = (2*pi*140./wavelengths).';
figure,
i=1;
hold on

clearvars max_energy leg
idx_ff = 1;

for idx_cr = 1:L1

for counter = 1:max_order*(max_order+1)
    [idxn, idxm] =ind2sub([max_order, max_order+1], counter);
    max_energy = Astats.energy.max(idx_cr, idx_ff, :, :, idxn, idxm);
    if max_energy > 1e-1
        %plot(wavelengths, (squeeze(Astats.mean.mag(:, :,idxn, idxm))))
        x = (2*pi*140./wavelengths).';
        plot(wavelengths, 2.*(2*idxn+1).*(squeeze(Astats.mean.mag(idx_cr, idx_ff,:,:,idxn, idxm)))./(x.^2))
        leg{i} = strjoin({'A n = ',num2str(idxn), ' m = ', num2str(idxm)}); 
        i = i+1;
        plot(wavelengths, 2.*(2*idxn+1).*(squeeze(Bstats.mean.mag(idx_cr, idx_ff,:,:,idxn, idxm)))./(x.^2))
        leg{i} = strjoin({'B n = ',num2str(idxn), ' m = ', num2str(idxm)}); 
        i = i+1;
    end
end
end
hold off  
legend(leg); 
    
%%
    
figure, 
i = 1;
hold on 
for counter = 1:max_order*(max_order+1)
    [idxn, idxm] =ind2sub([max_order, max_order+1], counter);
    max_energy = Astats.energy.max(idx_cr, idx_ff, :, :, idxn, idxm);
    if max_energy > 1e-1
        plot(wavelengths, (squeeze(Astats.mean.phase(idx_cr, idx_ff, :, :,idxn, idxm))))
        leg{i} = strjoin({'n = ',num2str(idxn), ' m = ', num2str(idxm)}); 
        i = i+1;
    end
end
hold off  
legend(leg); 

%end
    



%%
figure, 
hold on
for idx_r = 1:size(A,1)
plot(wavelengths, log10(squeeze(mean(Ipar(idx_r, 1, :, :,1), 3)./mean(Ipar(idx_r, 1, :, :,end), 3))));
end
hold off
%{

figure, 
plot(wavelengths, log10(Iper(:,1)./Iper(:,end)))

I = Ipar;
I(:,:,2) = Iper;
I = mean(I,3);

figure, 
plot(wavelengths, log10(I(:,1)./I(:,end)))
%}