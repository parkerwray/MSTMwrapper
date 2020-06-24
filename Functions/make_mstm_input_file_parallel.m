
function make_mstm_input_file_parallel(directory, fname, Nspheres, m_medium, k, beam_type, beam_width, alpha, beta, pol, near_field_cords, near_field_resolution)
inputfilename = fname;
outputfilename = strcat(inputfilename,'_output');
fid = fopen(strcat(directory,'/',inputfilename,'.inp'),'w');

fprintf(fid,'begin_comment\n');

fprintf(fid, strcat(directory,'\n'));

fprintf(fid,'end_comment\n');

%% Number of parallel processors per simulation
fprintf(fid, 'sm_number_processors\n');
fprintf(fid, '10\n');

%% FLAGS

fprintf(fid, 'write_sphere_data\n');
fprintf(fid, '1\n');

fprintf(fid, 'track_iterations\n');
fprintf(fid, '1\n');

fprintf(fid,'append_output_file\n');
fprintf(fid,'0\n');

fprintf(fid, 'run_print_file\n');
fprintf(fid, '\n');

fprintf(fid, 'calculate_scattering_coefficients\n');
fprintf(fid, '1\n');

fprintf(fid, 'calculate_t_matrix\n');
fprintf(fid, '1\n');

fprintf(fid, 'store_translation_matrix\n');
fprintf(fid, '1\n'); %% Setting this to 0 will have memory improvements but slower.

fprintf(fid, 'delta_scattering_angle_deg\n');
fprintf(fid, '1\n');

fprintf(fid, 'incident_or_target_frame\n');
fprintf(fid, '0\n');

fprintf(fid, 'normalize_scattering_matrix\n');
fprintf(fid, '0\n');

fprintf(fid, 'fixed_or_random_orientation\n');
fprintf(fid, '0\n');

fprintf(fid, 'azimuth_average_scattering_matrix\n');
fprintf(fid, '0\n');

fprintf(fid, 'calculate_near_field\n');
fprintf(fid, '0\n');

%% Solution Error Tollerance

fprintf(fid, 'mie_epsilon\n');
fprintf(fid, '1.0d-10\n');

fprintf(fid, 'translation_epsilon\n');
fprintf(fid, '1.0d-9\n');

fprintf(fid, 'solution_epsilon\n');
fprintf(fid, '1.0d-10\n');

fprintf(fid, 'max_number_iterations\n');
fprintf(fid, '100000\n');

fprintf(fid, 'iterations_per_correction\n');
fprintf(fid, '20\n');

%% File Location Specifications
fprintf(fid,'output_file\n');
fprintf(fid,strcat(outputfilename,'.dat\n'));

fprintf(fid, 'scattering_angle_file\n');
fprintf(fid, strcat(inputfilename,'_scat_angles.dat\n'));

%% Stuff related to physical position of the sphers
fprintf(fid,'number_spheres\n');
fprintf(fid,'%d\n',Nspheres);

fprintf(fid,'sphere_position_file\n'); % specify sphere position in file
fprintf(fid,strcat(inputfilename,'_sphere_file.pos\n'));

fprintf(fid,'length_scale_factor\n'); % wavenumber which scales all lengths to be unitless.
fprintf(fid, strcat(strrep(num2str(k,'%e'),'e','d'),'\n'));


%% Medium Refractive index and Cirality things
fprintf(fid, 'real_ref_index_scale_factor\n');
fprintf(fid, '1.0d0\n');
fprintf(fid, 'imag_ref_index_scale_factor\n');
fprintf(fid, '1.0d0\n');
fprintf(fid, 'real_chiral_factor\n');
fprintf(fid, '0.0d0\n');
fprintf(fid, 'imag_chiral_factor\n');
fprintf(fid, '0.0d0\n');
fprintf(fid, 'medium_real_ref_index\n');
fprintf(fid, strcat(strrep(num2str(real(m_medium),'%e'),'e','d'),'\n'));
fprintf(fid, 'medium_imag_ref_index\n');
fprintf(fid, strcat(strrep(num2str(imag(m_medium),'%e'),'e','d'),'\n'));
fprintf(fid, 'medium_real_chiral_factor\n');
fprintf(fid, '0.0d0\n');
fprintf(fid, 'medium_imag_chiral_factor\n');
fprintf(fid, '0.0d0\n');

%%
fprintf(fid, 'target_euler_angles_deg\n');
fprintf(fid, '0.0d0,0.0d0,0.0d0\n');






fprintf(fid, 'number_scattering_angles\n');
fprintf(fid, '360\n');



fprintf(fid, 'gaussian_beam_constant\n');
if beam_type == 1
    fprintf(fid, strcat(strrep(num2str(1./(beam_width.*k),'%e'),'e','d'),'\n'));
    if (1./(beam_width.*k))>0.2
        'Beam focus is too strong!'
    end
else 
    fprintf(fid, '0\n');
end


fprintf(fid, 'gaussian_beam_focal_point\n');
fprintf(fid, '0.0d0,0.0d0,0.0d0\n');



fprintf(fid, 'incident_azimuth_angle_deg\n');
fprintf(fid, strcat(strrep(num2str(alpha,'%e'),'e','d'),'\n'));

fprintf(fid, 'incident_polar_angle_deg\n');
fprintf(fid, strcat(strrep(num2str(beta,'%e'),'e','d'),'\n'));

% Scattering Stuff

fprintf(fid, 'scattering_coefficient_file\n');
fprintf(fid, strcat([inputfilename,'_scat_coeffs.dat\n']));

% Near Field Stuff %


%{
Near field values are calculated in a rectangular grid lying in the plane denoted
y−z plane -> 1 
z-x plane -> 2
x-y plane -> 3
%}
fprintf(fid, 'near_field_plane_coord\n');
fprintf(fid, '2\n');

fprintf(fid, 'near_field_plane_position\n');
fprintf(fid, '0.0d0\n');

% Spacing is scaled by k
fprintf(fid, 'near_field_plane_vertices\n');

fprintf(fid, num2str(floor(near_field_cords(1).*k),'%d'));
fprintf(fid, ' ');
fprintf(fid, num2str(floor(near_field_cords(2).*k),'%d'));
fprintf(fid, ' ');
fprintf(fid, num2str(ceil(near_field_cords(3).*k),'%d'));
fprintf(fid, ' ');
fprintf(fid, num2str(ceil(near_field_cords(4).*k),'%d'));
fprintf(fid, '\n');

fprintf(fid, 'spacial_step_size\n');
fprintf(fid, strcat(strrep(num2str(near_field_resolution,'%e'),'e','d'),'\n'));

fprintf(fid, 'polarization_angle_deg\n');
fprintf(fid, strcat(strrep(num2str(pol,'%e'),'e','d'),'\n'));

fprintf(fid, 'near_field_output_file\n');
fprintf(fid, strcat([inputfilename,'_nf.dat\n']));
fprintf(fid, 'near_field_output_data\n');
fprintf(fid, '2\n');
fprintf(fid, 'plane_wave_epsilon\n');
fprintf(fid, '1.0d-5\n');

fprintf(fid, 't_matrix_file\n');
fprintf(fid, strcat([inputfilename,'_tm.dat\n']));
fprintf(fid, 't_matrix_convergence_epsilon\n');
fprintf(fid, '1.d-9\n');



%% Set Sweep Files
fprintf(fid, 'new_run\n')

fprintf(fid,'sphere_position_file\n'); % specify sphere position in file
fprintf(fid,strcat(inputfilename,'_sphere_file.pos\n'));

fprintf(fid,'length_scale_factor\n'); % wavenumber which scales all lengths to be unitless.
fprintf(fid, strcat(strrep(num2str(k,'%e'),'e','d'),'\n'));

fprintf(fid,'output_file\n');
fprintf(fid,strcat(outputfilename,'.dat\n'));

fprintf(fid, 'scattering_coefficient_file\n');
fprintf(fid, strcat([inputfilename,'_scat_coeffs.dat\n']));








fprintf(fid, 'end_of_options\n');

fclose(fid);
end









