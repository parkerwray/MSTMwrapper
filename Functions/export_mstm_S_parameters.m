clearvars -except File FileName PathName
clc
% [File, FileName, PathName] = get_file();
% 
fid = fopen(File,'r');

% Skip useless info
for i = 1:7
    fgetl(fid);
end

% Grab number of spheres and total volume size parameter 
a = fscanf(fid, '%f');
num_spheres = a(1);
vol_size_param = a(2);

% Grab file that contrains sphere position info 
fgetl(fid);
sphere_position_file = fgetl(fid);

% Skip useless info
for i = 1:8
    fgetl(fid);
end

% File that lists the scattering angles that were calculated
file_recorded_scattering_angles = fgetl(fid);
fgetl(fid);
fgetl(fid);
Nscat_angles = fscanf(fid, '%f');

% Skip useless info
for i = 1:17
    fgetl(fid);
end

% File that lists the scattering coefficients
file_recorded_scattering_coeffs = fgetl(fid);

% File that lists the near field data and near field info
file_recorded_near_fields = fgetl(fid);
fgetl(fid);
near_field_type = fscanf(fid, '%f');
fgetl(fid);
a = fscanf(fid, '%f');
near_field_plane = a(1);
near_field_position = a(2);
fgetl(fid);
near_field_vertices = fscanf(fid, '%f');
near_field_spatial_resolution = fscanf(fid, '%f');

% Grab input polarization state
Input_polarization = fscanf(fid, '%f');

% Skip useless info
for i = 1:5
    fgetl(fid);
end

% Get sphere size parameter, location, index, and cross sections for
% particular lambda

a = fscanf(fid, '%f');



fclose(fid);



















