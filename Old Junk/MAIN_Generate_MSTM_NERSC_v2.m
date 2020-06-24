
%parentdir = uigetdir('$ME','Specify Save Directory');
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
r = 40;
ff = 0.3;
bounds = [100, 100, 50]./r;
dimension = 3;
giggles = 1;

[cord, bounds, a, am, Nspheres] = ...
    make_random_fcc_v2(r, ff, bounds, giggles, dimension);
lower_bound = bounds(1,:);
upper_bound = bounds(2,:);
make_spheres_in_rectangle(cord, r, lower_bound, upper_bound)

%% Get material refractive index for spheres
mat_file = 'Z:\Project - Dusty Plasma MURI\Material Data\Refractive Index Info\Ag_Yang.csv';
name = 'mat';
mat = pw_read_refinfo_mat(name, mat_file, 1);

%% Compile sphere [radius, x, y, z, n, k] into matrix 
clearvars spheres
wavelengths = 300:10:400;
k = 2*pi./wavelengths;
dummy = [r.*ones(size(cord,1),1), cord];
mat = index_in_range(mat, wavelengths, 1);

for idx = 1:length(wavelengths)
  spheres(:,:,idx) = [dummy,...
        mat.n(idx).*ones(size(cord,1),1),...
        mat.k(idx).*ones(size(cord,1),1)];
end
%%
medium_m = 1;
% Nspheres = size(spheres,1);
for idx = 1:length(wavelengths)
fname= strcat('mstm',... %GIVE A FILE NAME!
        '_lda_', num2str(wavelengths(idx)),'nm');
make_mstm_NERSC_job(parentdir, fname, spheres(:,:,idx), Nspheres, medium_m, k(idx),...
    beam_type, beam_waist,alpha, beta, pol, near_field_cords,...
    near_field_resolution);
end
  