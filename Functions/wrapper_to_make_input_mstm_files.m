
function wrapper_to_make_input_mstm_files(dir, fname, sphere,...
    Nspheres, medium_m, k, beam_type, beam_waist, alpha, beta,...
    pol, near_field_cords, near_field_resolution)
            
make_mstm_input_file_v2(dir, fname, Nspheres, medium_m, k, beam_type,...
    beam_waist, alpha, beta, pol, near_field_cords, near_field_resolution);
make_mstm_sphere_file(dir, fname, sphere);
make_mstm_scat_angle_file(dir, fname);

end


























