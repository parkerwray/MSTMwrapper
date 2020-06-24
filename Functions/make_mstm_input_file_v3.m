
function make_mstm_input_file_v3(directory,...
                        fname,...
                        mstm_flags,...
                        convergence,...
                        input_beam,...
                        near_field,...
                        mstm_input_params)
inputfilename = fname;
outputfilename = strcat(inputfilename,'_output');
fid = fopen(strcat(directory,'/',inputfilename,'.inp'),'w');


k = mstm_input_params.k;


fprintf(fid,'begin_comment\n');

fprintf(fid, strcat(directory,'\n'));

fprintf(fid,'end_comment\n');


%% Stuff related to physical position of the sphers
fprintf(fid,'number_spheres\n');
fprintf(fid,'%d\n',mstm_input_params.Nspheres);

fprintf(fid,'sphere_position_file\n'); % specify sphere position in file
fprintf(fid,strcat(inputfilename,'_sphere_file.pos\n'));

fprintf(fid,'length_scale_factor\n'); % wavenumber which scales all lengths to be unitless.
fprintf(fid, strcat(strrep(num2str(mstm_input_params.k,'%e'),'e','d'),'\n'));

fprintf(fid, 'write_sphere_data\n');
fprintf(fid, strcat(num2str(mstm_flags.write_sphere_data,'%d'),'\n'));

%% Output file things
fprintf(fid,'output_file\n');
fprintf(fid,strcat(outputfilename,'.dat\n'));
fprintf(fid,'append_output_file\n');
fprintf(fid,'0\n');
fprintf(fid, 'run_print_file\n');
fprintf(fid, '\n');


% Refractive index things
fprintf(fid, 'real_ref_index_scale_factor\n');
fprintf(fid, strcat(strrep(...
    num2str(mstm_input_params.real_ref_index_scale_factor,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'imag_ref_index_scale_factor\n');
fprintf(fid, strcat(strrep(...
    num2str(mstm_input_params.imag_ref_index_scale_factor,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'real_chiral_factor\n');
fprintf(fid, strcat(strrep(...
    num2str(mstm_input_params.real_chiral_factor,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'imag_chiral_factor\n');
fprintf(fid, strcat(strrep(...
    num2str(mstm_input_params.imag_chiral_factor,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'medium_real_ref_index\n');
fprintf(fid, strcat(strrep(...
    num2str(mstm_input_params.medium_real_ref_index,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'medium_imag_ref_index\n');
fprintf(fid, strcat(strrep(...
    num2str(mstm_input_params.medium_imag_ref_index,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'medium_real_chiral_factor\n');
fprintf(fid, strcat(strrep(...
    num2str(mstm_input_params.medium_real_chiral_factor,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'medium_imag_chiral_factor\n');
fprintf(fid, strcat(strrep(...
    num2str(mstm_input_params.medium_imag_chiral_factor,...
    '%e'),'e','d'),'\n'));

%%
fprintf(fid, 'target_euler_angles_deg\n');
fprintf(fid, '0.0d0,0.0d0,0.0d0\n');

% This value determines the order nmax for Mie theory of a single sphere.
% Setting this value to a negative interger (e.g., -#) forces the Mie
% theory code to only account for up to # orders.
fprintf(fid, 'mie_epsilon\n');
fprintf(fid, strcat(strrep(...
    num2str(convergence.mie_epsilon,...
    '%e'),'e','d'),'\n'));

% This value of epsilon determines the nmax for the T matrix. 
% You can play with this to play with nmax. 
fprintf(fid, 'translation_epsilon\n');
fprintf(fid, strcat(strrep(...
    num2str(convergence.translation_epsilon,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'solution_epsilon\n');
fprintf(fid, strcat(strrep(...
    num2str(convergence.solution_epsilon,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'max_number_iterations\n');
fprintf(fid, '100000\n');

fprintf(fid, 'store_translation_matrix\n');
fprintf(fid, strcat(num2str(mstm_flags.store_translation_matrix,'%d')...
    ,'\n')); %% Setting this to 0 will have memory improvements but slower.

fprintf(fid, 'sm_number_processors\n');
fprintf(fid, strcat(num2str(convergence.sm_number_processors,'%d')...
    ,'\n')); 

% fprintf(fid, 'near_field_translation_distance\n');
% fprintf(fid, '-1\n');
fprintf(fid, 'iterations_per_correction\n');
fprintf(fid, strcat(num2str(convergence.iterations_per_correction,'%d')...
    ,'\n')); 

fprintf(fid, 'scattering_angle_file\n');
fprintf(fid, strcat(inputfilename,'_scat_angles.dat\n'));
%{
Filename for a text file containing pairs of scattering directions, listed as one
pair per line per the format
                            ?1, ?1
                            ?2, ?2
                            ?NA , ?NA
The angles ?i and ?i appearing in the file are in degrees, with 0 ? ? ? 180 and 0 ? ? ? 360. The
delimeter can be space(s), a comma, or a tab. If the file contains fewer than NA pairs, the value of NA
is changed to the number of read pairs If NA is not specified in the input file the code will read all pairs
in the scattering angle file, and will set NA accordingly. The scattering_angle_file option is used
only for fixed orientation calculations, and angles for random orientation will be set solely through
selection of N? or ??. There is no default scattering angle file.
%}

fprintf(fid, 'number_scattering_angles\n');
fprintf(fid, '360\n');
fprintf(fid, 'delta_scattering_angle_deg\n');
fprintf(fid, '1\n');
fprintf(fid, 'incident_or_target_frame\n');
fprintf(fid, '0\n');
%{
Integer switch, relevant only for fixed orientation calculations. This selects
the reference coordinate frame upon which the scattering angles are based. When = 0, the frame
is based on the incident field so that ? = 0 is the forward scattering direction (i.e., the direction of
incident field propagation) and ? = ?, ? = 180 is the target z axis direction. When = 1, the frame
is based on the target frame, so that ? = 0 corresponds to the target z axis and ? = ?, ? = ? is
the incident field direction. This option affects only the association of the angles ? and ?, as listed in
the output along with the corresponding scattering matrix values, with a coordinate frame. The value
of the matrix for a set direction in fixed space will not be altered by the selection of the coordinate
frame: the matrix elements are defined using the standard convention under which the transverse
vector components of the scattered electric field are established using the plane formed by the incident
and scattered directions as a reference.
%}
fprintf(fid, 'normalize_scattering_matrix\n');
fprintf(fid, strcat(num2str(mstm_flags.normalize_scattering_matrix,'%d')...
    ,'\n'));

fprintf(fid, 'gaussian_beam_constant\n');
if input_beam.beam_type == 1
    fprintf(fid, strcat(strrep(num2str(1./(input_beam.gaussian_beam_width.*k),'%e'),'e','d'),'\n'));
    if (1./(input_beam.gaussian_beam_width.*mstm_input_params.k))>0.2
        'Beam focus is too strong!'
    end
else 
    fprintf(fid, '0\n');
end


fprintf(fid, 'gaussian_beam_focal_point\n');
fprintf(fid, '0.0d0,0.0d0,0.0d0\n');

fprintf(fid, 'fixed_or_random_orientation\n');
fprintf(fid, strcat(num2str(mstm_flags.fixed_or_random_orientation,'%d')...
    ,'\n'));

fprintf(fid, 'azimuth_average_scattering_matrix\n');
fprintf(fid, strcat(num2str(mstm_flags.azimuth_average_scattering_matrix,'%d')...
    ,'\n'));

fprintf(fid, 'incident_azimuth_angle_deg\n');
fprintf(fid, strcat(strrep(num2str(input_beam.incident_azimuth_angle_deg,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'incident_polar_angle_deg\n');
fprintf(fid, strcat(strrep(num2str(input_beam.incident_polar_angle_deg,...
    '%e'),'e','d'),'\n'));


% Scattering Stuff
fprintf(fid, 'calculate_scattering_coefficients\n');
fprintf(fid, strcat(num2str(mstm_flags.calculate_scattering_coefficients,'%d')...
    ,'\n'));

fprintf(fid, 'scattering_coefficient_file\n');
fprintf(fid, strcat([inputfilename,'_scat_coeffs.dat\n']));

% Near Field Stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid, 'track_iterations\n');
fprintf(fid, strcat(num2str(mstm_flags.track_nearfield_iterations,'%d')...
    ,'\n'));
fprintf(fid, 'calculate_near_field\n');
fprintf(fid, strcat(num2str(mstm_flags.calculate_near_field,'%d')...
    ,'\n'));

%{
Near field values are calculated in a rectangular grid lying in the plane denoted
y−z plane -> 1 
z-x plane -> 2
x-y plane -> 3
%}
fprintf(fid, 'near_field_plane_coord\n');
fprintf(fid, strcat(num2str(near_field.plane_cord,'%d')...
    ,'\n'));

fprintf(fid, 'near_field_plane_position\n');
fprintf(fid, strcat(strrep(num2str(near_field.plane_position,...
    '%e'),'e','d'),'\n'));

% Spacing is scaled by k
fprintf(fid, 'near_field_plane_vertices\n');

fprintf(fid, num2str(floor(near_field.plane_vertices(1).*k),'%d'));
fprintf(fid, ' ');
fprintf(fid, num2str(floor(near_field.plane_vertices(2).*k),'%d'));
fprintf(fid, ' ');
fprintf(fid, num2str(ceil(near_field.plane_vertices(3).*k),'%d'));
fprintf(fid, ' ');
fprintf(fid, num2str(ceil(near_field.plane_vertices(4).*k),'%d'));
fprintf(fid, '\n');


fprintf(fid, 'spacial_step_size\n');
fprintf(fid, strcat(strrep(num2str(near_field.resolution,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'polarization_angle_deg\n');
fprintf(fid, strcat(strrep(num2str(input_beam.polarization_angle_deg,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'near_field_output_file\n');
fprintf(fid, strcat([inputfilename,'_nf.dat\n']));

fprintf(fid, 'near_field_output_data\n');
fprintf(fid, '2\n');

fprintf(fid, 'plane_wave_epsilon\n');
fprintf(fid, strcat(strrep(num2str(convergence.plane_wave_epsilon,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'calculate_t_matrix\n');
fprintf(fid, strcat(num2str(mstm_flags.calculate_t_matrix,'%d')...
    ,'\n'));

fprintf(fid, 't_matrix_file\n');
fprintf(fid, strcat([inputfilename,'_tm.dat\n']));

fprintf(fid, 't_matrix_convergence_epsilon\n');
fprintf(fid, strcat(strrep(num2str(convergence.t_matrix_convergence_epsilon,...
    '%e'),'e','d'),'\n'));

fprintf(fid, 'end_of_options\n');

fclose(fid);
end









