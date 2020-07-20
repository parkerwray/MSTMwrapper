%% This code is used to process the results of 2D particle film simulations of random particle size distributions and random particle placement. The code extracts the modes of the center particle as a funciton of center particle size, cluster fill fraction, and wavelength. Statistics are generated beased on different unique distributions and polarization. The resulting statistics of the mode profile is then compared to the isolated particle case. 
% % 
% %% SPECIFY RELEVANT DIRECTORIES 

% Paths to import programs.
addpath(genpath('/home/elipaul/hypnos/codes/MSTMwrapper')); 
addpath(genpath('/home/elipaul/hypnos/codes/Matlab Functions'));  
addpath(genpath('/home/elipaul/hypnos/codes/randomparticles'));  

% Location of the fortran compiled run program for MSTM code.
mstm_location = '/home/elipaul/hypnos/Codes/MSTMwrapper/mstm_parallel_ubuntu.out';

% Specify the directory where you will place input and output files for
% communicating with the MSTM fortran program. Save files will also go to
% this location. 
parentdir = uigetdir('/home/elipaul/hypnos',' Specify Save Directory'); % for Linux
%parentdir = uigetdir('$ME','Specify Save Directory'); % for NERSC

% Change directory to location where input/output files are saved
oldFolder = cd(parentdir);

%% PROCESS MODE COEFFICIENTS 


clearvars A B Ax Bx size_param sphere_coffs_par sphere_coeffs_per Ipar Iper
max_order = 10;
theta = [0, pi];

parfor counter = 1:L1*L2*L3*L4*L5
    if mod(counter, 1000) == 0
        disp(counter);
        
    end
    [idx5, idx1, idx2, idx3, idx4] =ind2sub([L5, L1, L2, L3, L4], counter);
    fname =...
            strcat('mstm_',... 
            sprintf( '%03d', sub2ind([L1, L2, L3, L4, L5], idx1, idx2, idx3, idx4, idx5)));
        
    filename = strcat(parentdir,'/',fname);
    
    sphere_coeffs_data = cell(2, 1);
    % Convert linear index notation to vector index notation 
    for idx5 = 1:L5
        [sphere_coeffs_data{idx5}] =...
            export_mstm_scattering_coeffs_v3(strcat(filename,'_scat_coeffs.dat'));
    end

        sphere_coeffs_par = squeeze(sphere_coeffs_data(1));
        sphere_coeffs_per = squeeze(sphere_coeffs_data(2));
        
        [A(counter,:, :, :), B(counter,:, :, :),...
            Ax(counter, :, :, :), Bx(counter, :, :, :), size_param] = ...
                    process_convert_modes_v2(max_order, sphere_coeffs_par.', sphere_coeffs_per.');
          
        [Ipar(counter,:), Iper(counter,:)] = make_S_v2(A(counter,:, :, :), B(counter,:, :, :),...
            Ax(counter,:, :, :), Bx(counter,:, :, :), theta);  
  
end

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