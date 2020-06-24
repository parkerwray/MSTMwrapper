clear;
clc;

parentdir = uigetdir('$ME','Specify Save Directory');
% parentdir = uigetdir('Z:','Specify Save Directory');

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


%% Make Randomly positioned same sized spheres inside of a rectangle
clearvars r ff bounds dimension giggles cord a am Nspheres
clc
r = 82;
ff = 0.3;
bounds = [10000, 10000, 1];
%bounds = [1000, 1000, 1000];
dimension = 2;
%giggles = 500;
giggles = 5000;

[cord, bounds, a, am, Nspheres] = ...
    make_random_fcc_v3(r, ff, bounds, giggles, dimension);
lower_bound = bounds(1,:);
upper_bound = bounds(2,:);
make_spheres_in_rectangle(cord, r, lower_bound, upper_bound)

%% Get material refractive index for spheres
mat_file = '/global/cfs/cdirs/m2542/parkerwray/20_06_15_Material Data/Refractive Index Info/Ag_Yang.csv';
name = 'Ag';
mat1 = pw_read_refinfo_mat(name, mat_file, 1);

mat_file = '/global/cfs/cdirs/m2542/parkerwray/20_06_15_Material Data/Refractive Index Info/Ge.csv';
name = 'Ge';
mat2 = pw_read_refinfo_mat(name, mat_file, 1);
mat2.k(mat2.k==0) = 1.*10^(-6);

wavelengths = [400:50:500]; %[400:10:1000, 1100:100:5000, 5100:200:10000];
mat1 = index_in_range(mat1, wavelengths, 1);
mat2 = index_in_range(mat2, wavelengths, 1);
k = 2*pi./wavelengths;

%% Compile sphere [radius, x, y, z, n, k] into matrix for outermost shell
clearvars spheres spheres_layer1 spheres_layer2 spheres_layer3

r = 82;
dummy = [r.*ones(size(cord,1),1), cord];
for idx = 1:length(wavelengths)
  spheres_layer3(:,:,idx) = [dummy,...
        mat1.n(idx).*ones(size(cord,1),1),...
        mat1.k(idx).*ones(size(cord,1),1)];
end

r = 80;
dummy = [r.*ones(size(cord,1),1), cord];
for idx = 1:length(wavelengths)
  spheres_layer2(:,:,idx) = [dummy,...
        mat2.n(idx).*ones(size(cord,1),1),...
        mat2.k(idx).*ones(size(cord,1),1)];
end

r = 40;
dummy = [r.*ones(size(cord,1),1), cord];
for idx = 1:length(wavelengths)
  spheres_layer1(:,:,idx) = [dummy,...
        mat1.n(idx).*ones(size(cord,1),1),...
        mat1.k(idx).*ones(size(cord,1),1)];
end

spheres = [spheres_layer3; spheres_layer2; spheres_layer1];


%%
Nspheres = size(spheres,1); % Account for spheres inside spheres.
medium_m = 1;
for idx = 1:length(wavelengths)
fname= strcat('mstm',... %GIVE A FILE NAME!
        '_lda_', num2str(wavelengths(idx)),'nm');
make_mstm_NERSC_job(parentdir, fname, spheres(:,:,idx), Nspheres, medium_m, k(idx),...
    beam_type, beam_waist,alpha, beta, pol, near_field_cords,...
    near_field_resolution);
end

%%
nodes = '2';
time = '00:05:00';
jobs = idx;
make_mstm_SLURM_KNL_array_file(parentdir, nodes, time, jobs)



