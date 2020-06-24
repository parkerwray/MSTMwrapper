clc
clearvars -except sphere_file mat_file parentdir
%%
% USER DEFINED SAVE LOCATIONS
GET_SPHERE_DATA = 0;
GET_SPHERE_INDEX = 0;
GET_SAVE_DIRECTORY = 1;


if GET_SPHERE_DATA == 1
    addpath(genpath('/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/Parkers_diffusion_limited_aggregation_code')); %NP growth code uses user definced object classes. These classes need to be present in the search path to properly assign the class to the imported variables. 
    file_location = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/Parkers_diffusion_limited_aggregation_code';
    [sphere_file, FileName, PathName] = get_file(file_location);
    load(sphere_file);
end

if GET_SPHERE_INDEX == 1
    file_location = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/';
    [mat_file, FileName, PathName] = get_file(file_location);
%matfile = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/Refractive_Index_Info_Materials/Si_nk_200_800nm_Aspnes.csv';
end

if GET_SAVE_DIRECTORY == 1
    parentdir = uigetdir;
end




%% IMPORT MATERAIL DATA
%{
This section imports material data and defines the wavelengh(s) that will
be simulated. 
%}

name = 'mat';
mat = pw_read_refinfo_mat(name, mat_file, 1);
figure, 
hold on 
plot(mat.lda, real(mat.m));
plot(mat.lda, imag(mat.m));
hold off

%% SIMULATION INPUT PARAMETERS

beam_type = 0; % 0 = plane wave, 1 = Gaussian beam

if beam_type == 1
    beam_waist = ceil(1./(0.19.*k)); % Spot RADIUS at focus
    Zr = (pi.*beam_waist.^2)./lda; % Rayleigh range in [nm]
    disp(['The input beam diameter is ', num2str(2.*beam_waist), 'nm at focus'])
    disp(['The depth of field is ', num2str(2.*Zr), 'nm'])
    disp(['The unitless parameter is ', num2str(1./(k.*beam_waist)), ' (should be <= 0.2)'])
else
    beam_waist = 0;
end



%beam_width = lda;   % 0 = plane wave, 1/w*k = Gaussian beam, beam_width = width [nm];
alpha = 0; % input wave polar angle
beta = 0; % input wave venith angle
pol = 0; % polarization state
near_field_cords = [-3000,-3000,3000,6000]; % [-x ,-y ,x ,y] Coordinates [unit nm]
near_field_resolution = 0.5;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  GENERATE MSTM FILES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
This code generates the mstm input files that are used to run MSTM
%} 

cluster = {};
excitation = {};
sphere_coeffs = {};
sphere_result = {};
sphere = {};
%sphere = cell(length(lda),1);
wavelengths = 500:5:850;
sphere_cords = [0 0 0];
Nspheres = 1;
sphere_radi = 300;
medium_m = 1;

    %% Loop over wavelengths
for idx2 = 1:length(wavelengths)

    loud = 0;
    set_lda = wavelengths(idx2);
    k = 2*pi/(set_lda);
%     [sphere_m, medium_m, k, lda] = get_index(mat, set_lda, Nspheres, loud);


    sphere{idx2} = [sphere_radi, sphere_cords,...
        real(mat.m(idx2)),...
        imag(mat.m(idx2))];

    % Make directory to house the simulation files. Do this because the
    % Fortran code only accepts file names of a certain size, so you can be
    % more descriptive about the simulation from a folder!

    if beam_type == 1
        beam = ['GB_width_',num2str(beam_waist), 'nm_'];
    end
    if beam_type == 0
        beam = ['PW_'];
    end    

    % Calculate for both polarization states
    for idx3 = 1:2
        if idx3 == 1
            alpha = 0; %pol = 0; % polarization state
        else
            alpha = 90; %pol = 90;
        end
        dirname= strcat('Single_GaN_',... %GIVE A FILE NAME!
                'lda_', num2str(set_lda),'nm');

        mkdir(parentdir, dirname(1,:));
        dir = strcat(parentdir,'/',dirname(1,:));

        fname = 'mstm_sim';

        % Make mstm files to be used
        make_mstm_input_file_v2(dir, fname, Nspheres, medium_m, k, beam_type, beam_waist, alpha, beta, pol, near_field_cords, near_field_resolution);
        make_mstm_sphere_file(dir, fname, sphere{idx2});
        make_mstm_scat_angle_file(dir, fname);

        % Copy mstm files to mstm_fortran for simulation
        source_file = strcat(dir,'/',fname,'_single.inp');
        destination = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/mstm_fortran';
        copyfile(source_file, destination);

        source_file = strcat(dir, '/', fname,'_single_sphere_file.pos');
        destination = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/mstm_fortran';
        copyfile(source_file, destination);

        source_file = strcat(dir, '/', fname,'_single_scat_angles.dat');
        destination = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/mstm_fortran';
        copyfile(source_file, destination);

        % Run MSTM code
        command = './mstm_serial_ubuntu.exe mstm_sim_single.inp';
        %[status,cmdout] = system(command,'-echo');
        [status,cmdout] = system(command);

        % Load Results
        File_out = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/mstm_fortran/mstm_sim_output.dat';
        copyfile(File_out, dir);
        File_scat = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/mstm_fortran/mstm_sim_scat_coeffs.dat';
        copyfile(File_scat, dir);

        if alpha == 0 
            [cluster_par{idx2}, sphere_result_par{idx2}, excitation_par{idx2}] = export_output_file(File_out);
            [sphere_coeffs_par{idx2}] = export_mstm_scattering_coeffs(File_scat);
        elseif alpha == 90
            [cluster_perp{idx2}, sphere_result_perp{idx2}, excitation_perp{idx2}] = export_output_file(File_out);
            [sphere_coeffs_perp{idx2}] = export_mstm_scattering_coeffs(File_scat);
        end 
    end
    disp(['Wavelength ', num2str(wavelengths(idx2)), 'nm done! Time ', num2str(toc), 's']); 
end
disp(['Time to loop wavelengths ', num2str(toc/60),' min'])

    
    
    
save_file = strcat(parentdir,'/A_coeffs_results.mat');
save(save_file,...
    'sphere_coeffs_par', ...
    'sphere_coeffs_perp', ...
    'ff',...
    'ff_calc',...
    'wavelengths');



save_file = strcat(parentdir,'/A_full_results.mat');
save(save_file,...
    'sphere_coeffs_par', ...
    'sphere_result_par', ...
    'cluster_par',...
    'sphere_coeffs_perp', ...
    'sphere_result_perp', ...
    'cluster_perp',...
    'ff',...
    'ff_calc',...
    'I',...
    'wavelengths');


[modes_perp, modes_par, modes] = grab_modes(sphere_coeffs_perp, sphere_coeffs_par);
[A11_mag, B11_mag, Ax11_mag, Bx11_mag] = get_AB_mag_stats(modes);
[I0_par_c, I0_perp_c, I180_par_c, I180_perp_c] = get_I_compact(modes);
[qe,qs,qa, qds] = get_efficiencies(sphere_result_par);


save_file = strcat(parentdir,'/A_processed_results.mat');
save(save_file,...
    'modes_par', ...
    'modes_perp', ...
    'modes',...
    'A11_mag', ...
    'B11_mag', ...
    'Ax11_mag',...
    'Bx11_mag',...
    'I0_par_c',...
    'I0_perp_c',...
    'I180_par_c',...
    'I180_perp_c',...
    'qe',...
    'qs',...
    'qa',...
    'qds',...
    'wavelengths');

disp(['Saved!'])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  SUBFUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I2 = draw_spheres(I, cords, r)

for idx1 = 1:size(cords, 3)
    dummy = I;
    for idx2 = 1:size(cords, 1)
        circle = make_circle(I, cords(idx2,:,idx1)+round(([size(I,1), size(I,2),0]./2)), r);
        dummy = dummy+circle;
    end 
    I2(:,:,idx1) = dummy;
end

end

function circle = make_circle(I, cord, r)

imageSizeX = size(I,1);
imageSizeY = size(I,2);
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = cord(1);
centerY = cord(2);
radius = r;
circle= (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;

end

function [sphere_cords, Nspheres, sphere_radi] = get_random_sphere_positions(r,ff,sim_length, loud)

%% DEFINE SPHERE POSITION
% Spheres should be defined in the same unit as the wavelength. The MSTM
% code works with size parameter r/lda, which should be unitless.
% Distances between spheres should also be in the same unit system. 
% Sphere coordinates: (i,1) = x, (i,2) = y for the ith sphere

%   UPDATED CODE 10/20/2019
% The random sphere generation code has been updated to ensure fill
% fractions are correct. This code works for fill fractions up to 40%.
% Anything denser will not work. 

Asim = (sim_length)^2;
Asphere = pi*r^2; %% Old code was 2pi should be pi
Nspheres0 = round(ff*Asim/Asphere);
distance = (sim_length-2*r);
sphere_cords = [0,0,0];
idx = 1;
touching = 0;
no_room = 0;
tries = 0;
while idx < Nspheres0
    
    loc = [randi([-distance/2,distance/2],1,2),0];
    for idx2=1:size(sphere_cords,1)
        if norm(sphere_cords(idx2,:)-loc) < 2*r
            touching = 1;
        end
    end
    
    if touching == 0
        sphere_cords = [sphere_cords;loc];
        no_room = 0;
        idx = idx+1;
    end
    
    if touching == 1
        no_room = no_room+1;
    end
    
    if no_room >100
        disp('no more room!')
        sphere_cords = [];
        idx = 0;
        tries = tries+1;
        %break;
        
    end
    if tries == 100
        disp('not possible!')
        keyboard;
        %break;
    end
    
    touching = 0;
end
    
Nspheres = size(sphere_cords,1);
sphere_radi = r*ones(Nspheres,1); % implied nm
if loud == 1
    disp(['Number of spheres = ', num2str(Nspheres)])
    disp(['Fill fraction= ', num2str(Nspheres*Asphere/Asim)])
    disp(['Sphere radius distribution = ' num2str(mean(sphere_radi)), ' +/- ' num2str(std(sphere_radi)), 'nm'])

    x_min = min(sphere_cords(:,1))-sphere_radi(1);
    x_max = max(sphere_cords(:,1))+sphere_radi(1);
    disp(['X simulation size: ', num2str(x_min), ' x ' num2str(x_max) 'nm d = ' num2str(x_max-x_min), 'nm'])

    y_min = min(sphere_cords(:,2))-sphere_radi(1);
    y_max = max(sphere_cords(:,2))+sphere_radi(1);
    disp(['Y simulation size: ', num2str(y_min), ' x ' num2str(y_max) 'nm d = ' num2str(y_max-y_min), 'nm'])

    z_min = min(sphere_cords(:,3))-sphere_radi(1);
    z_max = max(sphere_cords(:,3))+sphere_radi(1);
    disp(['Z simulation size: ', num2str(z_min), ' x ' num2str(z_max) 'nm d = ' num2str(z_max-z_min), 'nm'])

end

end



function [sphere_m, medium_m, k, lda] = get_index(mat, set_lda, Nspheres, loud)

lda_min = ceil(min(mat.lda));
lda_max = floor(max(mat.lda));

lda = lda_min:1:lda_max;
n = interp1(mat.lda, mat.n, lda);
k = interp1(mat.lda, mat.k, lda);


idx = find(lda==set_lda);
lda = lda(idx); %EACH WAVELENGTH IS A NEW SIMULATION!
sphere_m = n(idx)+1i.*k(idx);  
sphere_m = repmat(sphere_m,Nspheres,1);
medium_m = 1+1i*0;
k = 2.*pi./lda;

if loud == 1
    figure, 
    plot(mat.lda, mat.n, mat.lda, mat.k)
    hold on 
    plot(lda, n, lda, k)
    xlabel('Wavelength (nm)')
end




end



function repeat = check_repeat(sphere_cords)

    last = length(sphere_cords);
    repeat = 0;
    for idx=1:last-1
        if isequal(sphere_cords{idx}, sphere_cords{last})
            repeat = 1;
            break;
        end
    end
    
end



function [A11_mag, B11_mag, Ax11_mag, Bx11_mag] = get_AB_mag_stats(modes)

A11_mag_avg = (squeeze(  mean(squeeze(abs(modes.Amn(2,:,:)).^2),2)  )); %./(modes.x(:,1).^2);
B11_mag_avg = (squeeze(  mean(squeeze(abs(modes.Bmn(2,:,:)).^2),2)  )); %./(modes.x(:,1).^2);
Ax11_mag_avg = (squeeze( mean(squeeze(abs(modes.Axmn(2,:,:)).^2),2) )); %./(modes.x(:,1).^2);
Bx11_mag_avg = (squeeze( mean(squeeze(abs(modes.Bxmn(2,:,:)).^2),2) )); %./(modes.x(:,1).^2);


A11_mag_std = (squeeze(  std(squeeze(abs(modes.Amn(2,:,:)).^2),[],2)  )); %./(modes.x(:,1).^2);
B11_mag_std = (squeeze(  std(squeeze(abs(modes.Bmn(2,:,:)).^2),[],2)  )); %./(modes.x(:,1).^2);
Ax11_mag_std = (squeeze( std(squeeze(abs(modes.Axmn(2,:,:)).^2),[],2) )); %./(modes.x(:,1).^2);
Bx11_mag_std = (squeeze( std(squeeze(abs(modes.Bxmn(2,:,:)).^2),[],2) )); %./(modes.x(:,1).^2);

A11_mag = [A11_mag_avg,A11_mag_std];
B11_mag = [B11_mag_avg,B11_mag_std];
Ax11_mag = [Ax11_mag_avg,Ax11_mag_std];
Bx11_mag = [Bx11_mag_avg,Bx11_mag_std];


end




function [qe,qs,qa, qds] = get_efficiencies(sphere_results)
    for lda_idx = 1:size(sphere_results,1)
        for instance_idx = 1:size(sphere_results,2)
            % Grab all the coeffs for a specific lda, instance pair
            qe(lda_idx, instance_idx) = sphere_results{lda_idx, instance_idx}.qe; 
            qs(lda_idx, instance_idx) = sphere_results{lda_idx, instance_idx}.qs;             
            qa(lda_idx, instance_idx) = sphere_results{lda_idx, instance_idx}.qa;      
        end
    end
qds = qe-qs-qa;


end

function [I0_par, I0_perp, I180_par, I180_perp] = get_I(modes_perp, modes_par)

    %% GET S AT THETA = 0, 180
    %[a011,a111,a021,a121,a221,a031,a131,a231,a331]
    a1n1_perp = squeeze(modes_perp.anm1([2,4,7],:,:));
    a1n2_perp = squeeze(modes_perp.anm2([2,4,7],:,:));
    b1n1_perp = squeeze(modes_perp.bnm1([2,4,7],:,:));
    b1n2_perp = squeeze(modes_perp.bnm2([2,4,7],:,:));

    a1n1_par = squeeze(modes_par.anm1([2,4,7],:,:));
    a1n2_par = squeeze(modes_par.anm2([2,4,7],:,:));
    b1n1_par = squeeze(modes_par.bnm1([2,4,7],:,:));
    b1n2_par = squeeze(modes_par.bnm2([2,4,7],:,:));
    %% GET S FOR THETA = 180 

    S1_180 = zeros(1,size(a1n1_perp,2));
    S2_180 = zeros(1,size(a1n1_perp,2));
    S3_180 = zeros(1,size(a1n1_perp,2));
    S4_180 = zeros(1,size(a1n1_perp,2));
    for n = 1:3
        S1_180 = S1_180+((-1i).^n).*((-1).^n).*(2*n+1).*(b1n1_perp(n,:,:)-a1n2_perp(n,:,:));
        S2_180 = S2_180+((-1i).^(n+1)).*((-1).^n).*(2*n+1).*(b1n2_par(n,:,:)-a1n1_par(n,:,:));
        S3_180 = S3_180+((-1i).^(n+1)).*((-1).^n).*(2*n+1).*(b1n2_perp(n,:,:)-a1n1_perp(n,:,:));
        S4_180 = S4_180+((-1i).^n).*((-1).^n).*(2*n+1).*(b1n1_par(n,:,:)-a1n2_par(n,:,:));
    end

    %% GET S FOR THETA = 0 

    S1_0 = zeros(1,size(a1n1_perp,2));
    S2_0 = zeros(1,size(a1n1_perp,2));
    S3_0 = zeros(1,size(a1n1_perp,2));
    S4_0 = zeros(1,size(a1n1_perp,2));
    for n = 1:3
        S1_0 = S1_0+((-1i).^n).*((-1)).*(2*n+1).*(b1n1_perp(n,:,:)+a1n2_perp(n,:,:));
        S2_0 = S2_0+((-1i).^(n+1)).*((-1)).*(2*n+1).*(b1n2_par(n,:,:)+a1n1_par(n,:,:));
        S3_0 = S3_0+((-1i).^(n+1)).*((-1)).*(2*n+1).*(b1n2_perp(n,:,:)+a1n1_perp(n,:,:));
        S4_0 = S4_0+((-1i).^n).*((-1)).*(2*n+1).*(b1n1_par(n,:,:)+a1n2_par(n,:,:));
    end
    %%
    S1_0 = squeeze(S1_0);
    S2_0 = squeeze(S2_0);
    S3_0 = squeeze(S3_0);
    S4_0 = squeeze(S4_0); 
    S1_180 = squeeze(S1_180);
    S2_180 = squeeze(S2_180);
    S3_180 = squeeze(S3_180);
    S4_180 = squeeze(S4_180);      
     
%     figure, 
%     plot(mean(abs(S1_0),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(abs(S2_0),2),':','LineWidth',6)   
%     plot(mean(abs(S3_0),2),'-','LineWidth',6)   
%     plot(mean(abs(S4_0),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
%     figure, 
%     plot(mean(abs(S1_180),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(abs(S2_180),2),':','LineWidth',6)   
%     plot(mean(abs(S3_180),2),'-','LineWidth',6)   
%     plot(mean(abs(S4_180),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
%     
%     figure, 
%     plot(mean(angle(S1_0),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(angle(S2_0),2),':','LineWidth',6)   
%     plot(mean(angle(S3_0),2),'-','LineWidth',6)   
%     plot(mean(angle(S4_0),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
    
    for idx = 1:size(S1_0,1)
        for idx2 = 1:size(S1_0,2)
        s23_iso(idx,idx2) = abs(S2_0(idx,idx2)).^2+abs(S3_0(idx,idx2)).^2;    
        s23_0_mix(idx,idx2) = dot(S2_0(idx,idx2),conj(S3_0(idx,idx2)))+dot(conj(S2_0(idx,idx2)),(S3_0(idx,idx2)));
        s23_180_mix(idx, idx2) = dot(S2_180(idx,idx2),conj(S3_180(idx,idx2)))+dot(conj(S2_180(idx,idx2)),(S3_180(idx,idx2)));
        s14_0_mix(idx,idx2) = dot(S4_0(idx,idx2),conj(S1_0(idx,idx2)))+dot(conj(S4_0(idx,idx2)),(S1_0(idx,idx2)));
        s14_180_mix(idx,idx2) =  dot(S4_180(idx,idx2),conj(S1_180(idx,idx2)))+dot(conj(S4_180(idx,idx2)),(S1_180(idx,idx2))); 
        
        I0_perp(idx,idx2) =  abs(S2_0(idx,idx2)).^2+abs(S3_0(idx,idx2)).^2+...
            dot(S2_0(idx,idx2),conj(S3_0(idx,idx2)))+dot(conj(S2_0(idx,idx2)),(S3_0(idx,idx2)));
        
        I0_par(idx,idx2) =  abs(S4_0(idx,idx2)).^2+abs(S1_0(idx,idx2)).^2+...
            dot(S4_0(idx,idx2),conj(S1_0(idx,idx2)))+dot(conj(S4_0(idx,idx2)),(S1_0(idx,idx2)));
        
        I180_perp(idx,idx2) =  abs(S2_180(idx,idx2)).^2+abs(S3_180(idx,idx2)).^2+...
            dot(S2_180(idx,idx2),conj(S3_180(idx,idx2)))+dot(conj(S2_180(idx,idx2)),(S3_180(idx,idx2)));
        
        I180_par(idx,idx2) =  abs(S4_180(idx,idx2)).^2+abs(S1_180(idx,idx2)).^2+...
            dot(S4_180(idx,idx2),conj(S1_180(idx,idx2)))+dot(conj(S4_180(idx,idx2)),(S1_180(idx,idx2)));    
        end
    end  
    
%     figure,
%     plot(mean(s23_0_mix,2))
%     hold on 
%     plot(mean(s23_180_mix,2))
%     plot(mean(s14_0_mix,2))
%     plot(mean(s14_180_mix,2)) 
%     hold off
%     legend('S_{23}(\theta = 0)','S_{23}(\theta = 180)','S_{14}(\theta = 0)','S_{14}(\theta = 180)')
   %keyboard;
   
    
    
end

function [I0_par, I0_perp, I180_par, I180_perp] = get_I_compact(modes)

    %% GET S AT THETA = 0, 180
    %[a011,a111,a021,a121,a221,a031,a131,a231,a331]
    Amn = squeeze(modes.Amn([2,4,7],:,:));
    Bmn = squeeze(modes.Bmn([2,4,7],:,:));
    Axmn = squeeze(modes.Axmn([2,4,7],:,:));
    Bxmn = squeeze(modes.Bxmn([2,4,7],:,:));

    %% GET S FOR THETA = 180 

    S1_180 = zeros(1,size(Amn,2));
    S2_180 = zeros(1,size(Amn,2));
    S3_180 = zeros(1,size(Amn,2));
    S4_180 = zeros(1,size(Amn,2));
    for n = 1:3
        S1_180 = S1_180+((-1i).^n).*((-1).^n).*(2*n+1).*(Bmn(n,:,:)-Amn(n,:,:));
        S2_180 = S2_180+((-1i).^(n+1)).*((-1).^n).*(2*n+1).*(Bmn(n,:,:)-Amn(n,:,:));
        S3_180 = S3_180+((-1i).^(n+1)).*((-1).^n).*(2*n+1).*(Bxmn(n,:,:)-Axmn(n,:,:));
        S4_180 = S4_180+((-1i).^n).*((-1).^n).*(2*n+1).*(Bxmn(n,:,:)-Axmn(n,:,:));
    end

    %% GET S FOR THETA = 0 

    S1_0 = zeros(1,size(Amn,2));
    S2_0 = zeros(1,size(Amn,2));
    S3_0 = zeros(1,size(Amn,2));
    S4_0 = zeros(1,size(Amn,2));
    for n = 1:3
        S1_0 = S1_0+((-1i).^n).*((-1)).*(2*n+1).*(Bmn(n,:,:)+Amn(n,:,:));
        S2_0 = S2_0+((-1i).^(n+1)).*((-1)).*(2*n+1).*(Bmn(n,:,:)+Amn(n,:,:));
        S3_0 = S3_0+((-1i).^(n+1)).*((-1)).*(2*n+1).*(Bxmn(n,:,:)+Axmn(n,:,:));
        S4_0 = S4_0+((-1i).^n).*((-1)).*(2*n+1).*(Bxmn(n,:,:)+Axmn(n,:,:));
    end
    %%
    S1_0 = squeeze(S1_0);
    S2_0 = squeeze(S2_0);
    S3_0 = squeeze(S3_0);
    S4_0 = squeeze(S4_0); 
    S1_180 = squeeze(S1_180);
    S2_180 = squeeze(S2_180);
    S3_180 = squeeze(S3_180);
    S4_180 = squeeze(S4_180);      
     
%     figure, 
%     plot(mean(abs(S1_0),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(abs(S2_0),2),':','LineWidth',6)   
%     plot(mean(abs(S3_0),2),'-','LineWidth',6)   
%     plot(mean(abs(S4_0),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
%     figure, 
%     plot(mean(abs(S1_180),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(abs(S2_180),2),':','LineWidth',6)   
%     plot(mean(abs(S3_180),2),'-','LineWidth',6)   
%     plot(mean(abs(S4_180),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
%     
%     figure, 
%     plot(mean(angle(S1_0),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(angle(S2_0),2),':','LineWidth',6)   
%     plot(mean(angle(S3_0),2),'-','LineWidth',6)   
%     plot(mean(angle(S4_0),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
    
    for idx = 1:size(S1_0,1)
        for idx2 = 1:size(S1_0,2)
        s23_iso(idx,idx2) = abs(S2_0(idx,idx2)).^2+abs(S3_0(idx,idx2)).^2;    
        s23_0_mix(idx,idx2) = dot(S2_0(idx,idx2),conj(S3_0(idx,idx2)))+dot(conj(S2_0(idx,idx2)),(S3_0(idx,idx2)));
        s23_180_mix(idx, idx2) = dot(S2_180(idx,idx2),conj(S3_180(idx,idx2)))+dot(conj(S2_180(idx,idx2)),(S3_180(idx,idx2)));
        s14_0_mix(idx,idx2) = dot(S4_0(idx,idx2),conj(S1_0(idx,idx2)))+dot(conj(S4_0(idx,idx2)),(S1_0(idx,idx2)));
        s14_180_mix(idx,idx2) =  dot(S4_180(idx,idx2),conj(S1_180(idx,idx2)))+dot(conj(S4_180(idx,idx2)),(S1_180(idx,idx2))); 
        
        I0_perp(idx,idx2) =  abs(S2_0(idx,idx2)).^2+abs(S3_0(idx,idx2)).^2+...
            dot(S2_0(idx,idx2),conj(S3_0(idx,idx2)))+dot(conj(S2_0(idx,idx2)),(S3_0(idx,idx2)));
        
        I0_par(idx,idx2) =  abs(S4_0(idx,idx2)).^2+abs(S1_0(idx,idx2)).^2+...
            dot(S4_0(idx,idx2),conj(S1_0(idx,idx2)))+dot(conj(S4_0(idx,idx2)),(S1_0(idx,idx2)));
        
        I180_perp(idx,idx2) =  abs(S2_180(idx,idx2)).^2+abs(S3_180(idx,idx2)).^2+...
            dot(S2_180(idx,idx2),conj(S3_180(idx,idx2)))+dot(conj(S2_180(idx,idx2)),(S3_180(idx,idx2)));
        
        I180_par(idx,idx2) =  abs(S4_180(idx,idx2)).^2+abs(S1_180(idx,idx2)).^2+...
            dot(S4_180(idx,idx2),conj(S1_180(idx,idx2)))+dot(conj(S4_180(idx,idx2)),(S1_180(idx,idx2)));    
        end
    end  
    
%     figure,
%     plot(mean(s23_0_mix,2))
%     hold on 
%     plot(mean(s23_180_mix,2))
%     plot(mean(s14_0_mix,2))
%     plot(mean(s14_180_mix,2)) 
%     hold off
%     legend('S_{23}(\theta = 0)','S_{23}(\theta = 180)','S_{14}(\theta = 0)','S_{14}(\theta = 180)')
%    keyboard;
   
    
    
end


function [modes_perp, modes_par, modes] = grab_modes(sphere_coeffs_perp, sphere_coeffs_par)

[anm1_par, anm2_par, bnm1_par, bnm2_par,x]= convert_modes(sphere_coeffs_par);
modes_par = struct(...
    'anm1', anm1_par,...
    'bnm1', bnm1_par,...
    'anm2', anm2_par,...
    'bnm2', bnm2_par,...
    'x',x);
 
[anm1_perp, anm2_perp, bnm1_perp, bnm2_perp,x]= convert_modes(sphere_coeffs_perp);
modes_perp = struct(...
    'anm1', anm1_perp,...
    'bnm1', bnm1_perp,...
    'anm2', anm2_perp,...
    'bnm2', bnm2_perp,...
    'x',x);

% modes = struct(...
%     'Amn', cat(3, anm2_perp, anm1_par),...
%     'Bmn', cat(3, bnm1_perp, bnm2_par),...
%     'Axmn', cat(3,anm1_perp, anm2_par),...
%     'Bxmn', cat(3,bnm2_perp, bnm1_par),...
%     'x',cat(2,x,x) );

modes = struct(...
    'Bmn', cat(3, anm2_perp, bnm2_par),...
    'Amn', cat(3, bnm1_perp, anm1_par),...
    'Bxmn', cat(3,anm1_perp, bnm1_par),...
    'Axmn', cat(3,bnm2_perp, anm2_par),...
    'x',cat(2,x,x) );



end

function [anm1, anm2, bnm1, bnm2, x,...
    qext_par, qsca_par, qabs_par,...
    qext_per, qsca_per, qabs_per]= convert_modes(sphere_coeffs)

%[-22,-12,02,12,22]
%[22,12,02,-12,-22]
%n = 2
% take [n+1,end] of both arrays
% add acording to m

dummy_te = {};
dummy_tm = {};

    for lda_idx = 1:size(sphere_coeffs,1)
        for instance_idx = 1:size(sphere_coeffs,2)
            % Grab all the coeffs for a specific lda, instance pair
            dummy_te = sphere_coeffs{lda_idx, instance_idx}.a_te;
            dummy_tm = sphere_coeffs{lda_idx, instance_idx}.a_tm;   
            x(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.ka;
            qsca_par(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qsca_par; 
            qabs_par(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qabs_par;             
            qext_par(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qext_par;             
            qsca_per(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qsca_per; 
            qabs_per(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qabs_per;             
            qext_per(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qext_per; 
            
            
            
            for n = 1:length(dummy_te)
                if n == 1
                   % [-11,01,11]
                   dummy_tm_n = dummy_te{n};
                   dummy_te_n = dummy_tm{n};

                   a111 = (dummy_tm_n(3)-dummy_tm_n(1));
                   a011 = dummy_tm_n(2);

                   a112 = 2.*(dummy_te_n(3)-dummy_te_n(1));
                   a012 = dummy_te_n(2); 

                   b111 = (dummy_tm_n(3)+dummy_tm_n(1));
                   b011 = dummy_tm_n(2);

                   b112 = 2.*(dummy_te_n(3)+dummy_te_n(1));
                   b012 = dummy_te_n(2);               
                end

                if n == 2
                   %[-22,-12,02,12,22]                
                   dummy_tm_n = dummy_te{n};
                   dummy_te_n = dummy_tm{n};

                   a221 = (dummy_tm_n(5)+dummy_tm_n(1));
                   a121 = (dummy_tm_n(4)-dummy_tm_n(2));              
                   a021 = dummy_tm_n(3);

                   a222 = 2.*(dummy_te_n(5)+dummy_te_n(1));
                   a122 = 2.*(dummy_te_n(4)-dummy_te_n(2));              
                   a022 = dummy_te_n(3);               


                   b221 = (dummy_tm_n(5)-dummy_tm_n(1));
                   b121 = (dummy_tm_n(4)+dummy_tm_n(2));              
                   b021 = dummy_tm_n(3);

                   b222 = 2.*(dummy_te_n(5)-dummy_te_n(1));
                   b122 = 2.*(dummy_te_n(4)+dummy_te_n(2));              
                   b022 = dummy_te_n(3);              
                end            

                if n == 3
                   %[-33,-23,-13,03,13,23,33]                
                   dummy_tm_n = dummy_te{n};
                   dummy_te_n = dummy_tm{n};

                   a331 = (dummy_tm_n(7)-dummy_tm_n(1));
                   a231 = (dummy_tm_n(6)+dummy_tm_n(2));
                   a131 = (dummy_tm_n(5)-dummy_tm_n(1));              
                   a031 = dummy_tm_n(4);

                   a332 = 2.*(dummy_te_n(7)-dummy_te_n(1));
                   a232 = 2.*(dummy_te_n(6)+dummy_te_n(2));
                   a132 = 2.*(dummy_te_n(5)-dummy_te_n(1));              
                   a032 = dummy_te_n(4);              

                   b331 = (dummy_tm_n(7)+dummy_tm_n(1));
                   b231 = (dummy_tm_n(6)-dummy_tm_n(2));
                   b131 = (dummy_tm_n(5)+dummy_tm_n(1));              
                   b031 = dummy_tm_n(4);

                   b332 = 2.*(dummy_te_n(7)+dummy_te_n(1));
                   b232 = 2.*(dummy_te_n(6)-dummy_te_n(2));
                   b132 = 2.*(dummy_te_n(5)+dummy_te_n(1));              
                   b032 = dummy_te_n(4);                             
                end               
            end

            anm1(:,lda_idx,instance_idx) = [a011,a111,a021,a121,a221,a031,a131,a231,a331]; 
            anm2(:,lda_idx,instance_idx) = [a012,a112,a022,a122,a222,a032,a132,a232,a332]; 
            bnm1(:,lda_idx,instance_idx) = [b011,b111,b021,b121,b221,b031,b131,b231,b331]; 
            bnm2(:,lda_idx,instance_idx) = [b012,b112,b022,b122,b222,b032,b132,b232,b332];         

        end
    end

end

%%%% This code makes dense packed spheres
function [cords2, Nspheres, sphere_radi] = dense_packed_cords(sim_length, r, ff)


[Nspheres, distance] = number_of_spheres(sim_length, r, ff);

cords = zeros(Nspheres,3);

cords(1,:) = -(distance/2).*[1,1,0];
flag = 0;
for idx=2:Nspheres
    
    next_loc = cords(idx-1,:)+[2*r+r/4, 0,0];
    dist = sqrt(next_loc(1).^2+next_loc(2).^2+next_loc(3).^3);
    
    if dist<2*r
        next_loc = [0,0,0]; %Make origin NP
        flag = 1;
    elseif next_loc(1)>distance/2
        cords(idx,:) = [-(distance/2), cords(idx-1,2)+2*r+r/4,0];
    else
        cords(idx,:) = next_loc;
    end

end

if flag == 0
    cords(Nspheres,:) = [0,0,0];
end

%visualize_spheres(cords, sim_length,r)


cords = flipud(cords);
for idx = 1:500
    for idx2 = 1:Nspheres

        me = cords(idx2,:);
        if me(1) == 0 && me(2) == 0 && me(3) == 0
            orig_idx = idx2;
            continue;
        end
        dummy_cords = cords;        
        dummy_cords(idx2,:) = [];
     
        new_me = giggle(me, dummy_cords, sim_length,r);
        cords(idx2,:) = new_me;
        

    end
end
cords2 = circshift(cords,Nspheres-orig_idx+1,1);
sphere_radi = r*ones(Nspheres,1); % implied nm
   
%visualize_spheres(cords2, sim_length,r)
end

function hit = check_me(me, cords, r, sim_length)

hit = 0; % flag goes 1 if you hit another sphere

% Check if sphere went over boundary
if abs(me(1)) > (sim_length-2*r)/2
    hit = 1;
elseif abs(me(2)) > (sim_length-2*r)/2
    hit = 1;
elseif abs(me(3)) > (sim_length-2*r)/2
    hit = 1;
end

% Check if sphere hit another sphere
if hit ~= 1

    for idx = 1:size(cords,1)

        dx = cords(idx, 1)-me(1);
        dy = cords(idx, 2)-me(2);
        dz = cords(idx, 3)-me(3);
        dr = sqrt(dx^2+dy^2+dz^2);

        if dr == 0 % Do nothing you are comparting to yourself
        elseif dr < 2.*r % you hit another sphere! Retry.
            hit = 1;
        end

        if hit == 1
            break
        end
    end

end

end

        
function me = giggle(me, cords, sim_length,r)

for counter = 1:100
    new_me = me + (r/10).*[randi([-1,1],1,2),0];
    hit = check_me(new_me, cords, r, sim_length);
    if hit == 0
        me = new_me;
    end
end


end



function [Nspheres, distance] = number_of_spheres(sim_length, r, ff)

% Get the number of spheres necessary to fill this required fill fraction

Asim = (sim_length)^2;
Asphere = pi*r^2; 
Nspheres = round(ff*Asim/Asphere);
distance = (sim_length-2*r);

end

% 
% function I2 = draw_spheres(I, cords, r)
% 
% for idx1 = 1:size(cords, 3)
%     dummy = I;
%     for idx2 = 1:size(cords, 1)
%         circle = make_circle(I, cords(idx2,:,idx1)+round(([size(I,1), size(I,2),0]./2)), r);
%         dummy = dummy+circle;
%     end 
%     I2(:,:,idx1) = dummy;
% end
% 
% end
% 
% function circle = make_circle(I, cord, r)
% 
% imageSizeX = size(I,1);
% imageSizeY = size(I,2);
% [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% % Next create the circle in the image.
% centerX = cord(1);
% centerY = cord(2);
% radius = r;
% circle= (rowsInImage - centerY).^2 ...
%     + (columnsInImage - centerX).^2 <= radius.^2;
% 
% end

function visualize_spheres(cords, sim_length,r)

    I = zeros(sim_length,sim_length);
    I = draw_spheres(I, cords, r);

    disp(['Fill Fraction ',num2str(100.*sum(sum(mean(I,3)))./(size(I,1).*size(I,2))), '%'])

    figure, 
    imshow(I)


end































