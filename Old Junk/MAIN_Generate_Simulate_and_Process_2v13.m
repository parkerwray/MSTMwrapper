clc
clearvars -except sphere_file mat_file parentdir
%%
% USER DEFINED SAVE LOCATIONS
GET_SPHERE_DATA = 0;
GET_SPHERE_INDEX = 1;
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



%%
loud = 1;
r = 44;
ff = 0.1;
sim_length = 20*r;
%[sphere_cords, Nspheres, sphere_radi] = get_random_sphere_positions(r,ff,sim_length, loud);


%% IMPORT MATERAIL DATA
%{
This section imports material data and defines the wavelengh(s) that will
be simulated. 
%}

name = 'mat';
mat = pw_read_refinfo_mat(name, mat_file, 0);



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

% sphere = struct{...
%     'index', sphere_m,...
%     'cords', sphere_cords,...
%     'N', Nspheres};
% 
% simulation = struct{...
%     'Input_Beam', [],...
%     'Beam_Width', beam_width,...
%     'Polar Angle', alpha,...
%     'Zenith_Angle', beta,...
%     'Polarization', pol,...
%     'Env_Index', medium_m,...
%     'Sphere_Index', sphere_m,...
%     'nf_cords', [-3000,-3000,3000,6000],...   % Unit nm
%     'nf_res', 0.5};
    


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
wavelengths = 370:2:500;
length_itterations = 25;
for idx = 1:25
    tic
    % Define parameters for random sphere generation
    loud = 1;
    
    
    %% Generate random spheres with sphere at origin
    counter = 0;
    if idx > 1
        repeat = 1;
        % Check to make sure random orientation is new!
        while repeat
            [sphere_cords{idx}, Nspheres, sphere_radi] = dense_packed_cords(sim_length, r, ff);
            visualize_spheres(sphere_cords{idx}, sim_length,r);
            %[sphere_cords{idx}, Nspheres, sphere_radi] = get_random_sphere_positions(r,ff,sim_length, loud);
            repeat = check_repeat(sphere_cords); %If no new, retry
            counter = counter+repeat;
            if counter > 1000
                keyboard; %No more unique orientations!
            end
        end
    else
        [sphere_cords{idx}, Nspheres, sphere_radi] = dense_packed_cords(sim_length, r, ff);
%         visualize_spheres(sphere_cords{idx}, sim_length,r);
        %[sphere_cords{idx}, Nspheres, sphere_radi] = get_random_sphere_positions(r,ff,sim_length, loud); %first random has nothing to compare
    end
    
    disp(['Time to generate sphere ', num2str(toc/60),' min'])
    
    %% Visualize the spheres and ensure correlation/fill fractions are correct
    cords_dummy = sphere_cords{idx};
    cords_dummy = reshape(cords_dummy,Nspheres,3,1);

    I(:,:,idx) = zeros(sim_length,sim_length);
    I(:,:,idx) = draw_spheres(I(:,:,idx), cords_dummy, r);
    ff_calc(idx) = (sum(sum(I(:,:,idx)))./(size(I,1).*size(I,2)));
    
    disp(['Fill Fraction ',num2str(100.*ff_calc(idx)), '%'])
    if (ff_calc(idx)-ff)./ff > 0.1
        keyboard;
    end
    
    tic
    %% Loop over wavelengths
    for idx2 = 1:length(wavelengths)

        loud = 0;
        set_lda = wavelengths(idx2);
        [sphere_m, medium_m, k, lda] = get_index(mat, set_lda, Nspheres, loud);


        sphere{idx,idx2} = [sphere_radi, sphere_cords{idx},...
            real(sphere_m),...
            imag(sphere_m)];

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
            dirname= strcat('mstm_Si_',... %GIVE A FILE NAME!
                    'Size_', sprintf('%02d',round((sim_length)./sphere_radi(1))),'_',...
                    'FF_', sprintf('%02d',100*ff),'_',...
                    'NP_',sprintf('%02d',Nspheres),'_',...
                    'lda_', sprintf('%02d',lda),'nm_',...
                     beam,...
                    'alpha_', sprintf('%02d',alpha),'deg_',...
                    'beta_', sprintf('%02d',beta),'deg_',...
                    'polarization_', sprintf('%02d',pol),'_',...
                    'radius_',num2str(sphere_radi(1)),'nm',...
                    '_run_',sprintf('%02d',idx));

            mkdir(parentdir, dirname(1,:));
            dir = strcat(parentdir,'/',dirname(1,:));

            fname = 'mstm_sim';

            % Make mstm files to be used
            make_mstm_input_file_v2(dir, fname, Nspheres, medium_m, k, beam_type, beam_waist, alpha, beta, pol, near_field_cords, near_field_resolution);
            make_mstm_sphere_file(dir, fname, sphere{idx,idx2});
            make_mstm_scat_angle_file(dir, fname);

            % Copy mstm files to mstm_fortran for simulation
            source_file = strcat(dir,'/',fname,'.inp');
            destination = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/mstm_fortran';
            copyfile(source_file, destination);

            source_file = strcat(dir, '/', fname,'_sphere_file.pos');
            destination = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/mstm_fortran';
            copyfile(source_file, destination);

            source_file = strcat(dir, '/', fname,'_scat_angles.dat');
            destination = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/mstm_fortran';
            copyfile(source_file, destination);

            % Run MSTM code
            command = '/usr/lib64/openmpi/bin/mpirun -n 10 ./mstm_parallel_ubuntu.out mstm_sim.inp';
            %[status,cmdout] = system(command,'-echo');
            [status,cmdout] = system(command);

            % Load Results
            File_out = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/mstm_fortran/mstm_sim_output.dat';
            copyfile(File_out, dir);
            File_scat = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/mstm_fortran/mstm_sim_scat_coeffs.dat';
            copyfile(File_scat, dir);
            
            if alpha == 0 
                [cluster_par{idx2,idx}, sphere_result_par{idx2,idx}, excitation_par{idx2,idx}] = export_output_file(File_out);
                [sphere_coeffs_par{idx2,idx}] = export_mstm_scattering_coeffs_v3(File_scat);
            elseif alpha == 90
                [cluster_perp{idx2,idx}, sphere_result_perp{idx2,idx}, excitation_perp{idx2,idx}] = export_output_file(File_out);
                [sphere_coeffs_perp{idx2,idx}] = export_mstm_scattering_coeffs_v3(File_scat);
            end 
        end
        disp(['Wavelength ', num2str(wavelengths(idx2)), 'nm done! Time ', num2str(toc), 's']); 
    end
    disp(['Time to loop wavelengths ', num2str(toc/60),' min'])
    
%     [A(idx,:,:,:,:), B(idx,:,:,:,:), Ax(idx,:,:,:,:), Bx(idx,:,:,:,:), size_param(idx,:)] =...
%         process_convert_modes(sphere_coeffs_par{:,idx}, sphere_coeffs_perp{:,idx});
    
end

save_file = strcat(parentdir,'/A_sim_data.mat');
save(save_file,...
    'I', ...
    'r',...
    'sim_length',...
    'mat', ...
    'ff',...
    'ff_calc',...
    'wavelengths');

save_file = strcat(parentdir,'/A_coeffs_results.mat');
save(save_file,...
    'sphere_coeffs_par', ...
    'sphere_coeffs_perp', ...
    'ff',...
    'ff_calc',...
    'wavelengths');

save_file = strcat(parentdir,'/A_sphere_results.mat');
save(save_file,...
    'sphere_result_par', ...
    'sphere_result_perp', ...
    'ff',...
    'ff_calc',...
    'wavelengths');

save_file = strcat(parentdir,'/A_cluster_results.mat');
save(save_file,...
    'cluster_par',...
    'cluster_perp',...
    'ff',...
    'ff_calc',...
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
for idx = 1:1000
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































