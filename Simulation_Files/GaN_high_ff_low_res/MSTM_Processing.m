%% This code is used to process the results of 2D particle film simulations of random particle size distributions and random particle placement. The code extracts the modes of the center particle as a funciton of center particle size, cluster fill fraction, and wavelength. Statistics are generated beased on different unique distributions and polarization. The resulting statistics of the mode profile is then compared to the isolated particle case. 
% 


%% PROCESS MODE COEFFICIENTS 



max_order = 10;
theta = [0, pi];


%%
clearvars A B Ax Bx size_param sphere_coffs_par sphere_coeffs_per Ipar Iper
for counter = 1:L1*L2*L3*L4
    
    % Convert linear index notation to vector index notation 
    [idx1, idx2, idx3, idx4] =ind2sub([L1, L2, L3, L4], counter);
    
    
    sphere_coeffs_par = squeeze(sphere_coeffs_data(idx1, idx2, idx3, idx4, 1));
    sphere_coeffs_per = squeeze(sphere_coeffs_data(idx1, idx2, idx3, idx4, 2));
    
    [A(idx1, idx2, idx3, idx4,:, :, :), B(idx1, idx2, idx3, idx4,:, :, :),...
        Ax(idx1, idx2, idx3, idx4, :, :, :), Bx(idx1, idx2, idx3, idx4, :, :, :), size_param] = ...
                process_convert_modes_v2(max_order, sphere_coeffs_par.', sphere_coeffs_per.');
      
    [Ipar(idx1, idx2, idx3, idx4,:), Iper(idx1, idx2, idx3, idx4,:)] = make_S_v2(A(idx1, idx2, idx3, idx4,:, :, :), B(idx1, idx2, idx3, idx4,:, :, :),...
        Ax(idx1, idx2, idx3, idx4,:, :, :), Bx(idx1, idx2, idx3, idx4,:, :, :), theta);  
      
    
    [qe(idx1, idx2, idx3, idx4), qa(idx1, idx2, idx3, idx4),...
        qsi(idx1, idx2, idx3, idx4), qsd(idx1, idx2, idx3, idx4)]=...
        get_efficiencies(sphere_coeffs_par);
    
    
    %{
        process_convert_modes was written under the assumption that the inputs
        were:
             1) The sphere coefficients for parallel and perpendicular
             illumination.
             2) Sphere coefficients are given as a {lda, Ndistributions} cell.
             3) All modes above max_order are negligable. 
      %}
  
end


%%
clearvars sphere_coeffs_data sphere_data cluster_data command fname sphere_coeffs_par sphere_coeffs_per
save('A_Modes_Far_Fields_Efficiencies.mat', 'A','B','Ax','Bx','Ipar','Iper','qe','qa','qsi','qsd');
%save('All_non_cell_data.mat')



%% GET MODE STATISTICS

clc
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
clc

figure,
i=1;
hold on

clearvars max_energy leg
idx_ff = 1;

for idx_cr = 1: L1 %L1

for counter = 1:max_order*(max_order+1)
    [idxn, idxm] =ind2sub([max_order, max_order+1], counter);
    max_energy = Astats.energy.max(idx_cr, idx_ff, :, :, idxn, idxm);
    if max_energy > 1e-1
        
        x = (2*pi*center_radiis(idx_cr)./wavelengths).';
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


    


%%
figure, 

I = squeeze((mean(Ipar,3)+mean(Iper,3))./2);
fbr = (I(:, :, :, 1)./I(:, :, :, end)); 



hold on
idx_ff = 1;
for idx_r = 1:5
plot(wavelengths, log10(squeeze(fbr(idx_r, idx_ff,:))) );
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