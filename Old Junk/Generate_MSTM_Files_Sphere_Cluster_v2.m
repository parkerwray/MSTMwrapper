% clear;
% clc;


%% IMPORT SPHERE DATA
%{ 
This secton of code is used to import sphere data from the nanoparticle 
deposition simulation 
%}
% addpath('Z:\Codes\diffusion-master\Parkers_diffusion_limited_aggregation_code'); %NP growth code uses user definced object classes. These classes need to be present in the search path to properly assign the class to the imported variables. 
% load('Z:\Project - Dusty Plasma MURI\19_02_18_Testing_Lumerical_Random_NP_Deposition_Simulations\Test_Large_NP_Simulaiton\Mid_Deposition_1957NPs.mat');

sphere_cords = space2.np_coordinates(:,:);
sphere_cords = sphere_cords - repmat(mean(sphere_cords,1), size(sphere_cords,1),1); % Normalize the coordinates such that the cluser is centered at zero.
Nspheres = size(sphere_cords,1);
sphere_radi = 44.*ones(Nspheres,1);

%% IMPORT MATERAIL DATA
%{
This section imports material data and defines the wavelengh(s) that will
be simulated. 
%}
matfile = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/Refractive_Index_Info_Materials/Si_nk_200_800nm_Aspnes.csv';
name = 'Si';
mat = pw_read_refinfo_mat(name, matfile,1);

lda_min = ceil(min(mat.lda));
lda_max = floor(max(mat.lda));

lda = lda_min:1:lda_max;
n = interp1(mat.lda, mat.n, lda);
k = interp1(mat.lda, mat.k, lda);

figure, 
plot(mat.lda, mat.n, mat.lda, mat.k)
hold on 
plot(lda, n, lda, k)
xlabel('Wavelength (nm)')

idx = find(lda==459);
lda = lda(idx); %backward = 32  forward = 34; %forward %EACH WAVELENGTH IS A NEW SIMULATION!
sphere_m = n(idx)+1i.*k(idx);  
sphere_m = repmat(sphere_m,Nspheres,1);

medium_m = 1+1i*0;

k = 2.*pi./lda;



%% GENERATE MSTM FILES
%{
This code generates the mstm input files that are used to run MSTM
%} 
parentdir = uigetdir;
sphere = cell(length(lda),1);
for idx = 1:length(lda)
sphere{idx} = [sphere_radi, sphere_cords,...
    real(sphere_m(:,idx)),...
    imag(sphere_m(:,idx))];

    
    [max(sphere_cords(:,2))   max(sphere_cords(:,2))].*k;

    % Make directory to house the simulation files. Do this because the
    % Fortran code only accepts file names of a certain size, so you can be
    % more descriptive about the simulation from a folder!
    dirname= strcat('mstm_2839NP_Si_Huygens_GB_Beta_20deg_r',... %GIVE A FILE NAME!
        num2str(sphere_radi(1)), 'nm_lda',...
        num2str(lda),'nm');
    mkdir(parentdir, dirname(1,:));
    dir = strcat(parentdir,'/',dirname(1,:));
    
    fname = 'mstm_2k_Si_F_H_GB_20deg';

    make_mstm_input_file(dir,fname,Nspheres,medium_m, k)
    make_mstm_sphere_file(dir,fname,sphere{idx})
    make_mstm_scat_angle_file(dir, fname)
end


















