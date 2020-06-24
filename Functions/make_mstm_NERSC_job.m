

function make_mstm_NERSC_job(dir, fname,sphere, Nspheres, medium_m, k, beam_type, beam_waist, alpha, beta, pol, near_field_cords, near_field_resolution)
            
% mkdir(parentdir, dirname(1,:));
% dir = strcat(parentdir,'/',dirname(1,:));
% fname = 'mstm_sim';
% Make mstm files to be used
make_mstm_input_file_v2(dir, fname, Nspheres, medium_m, k, beam_type,...
    beam_waist, alpha, beta, pol, near_field_cords, near_field_resolution);
make_mstm_sphere_file(dir, fname, sphere);
make_mstm_scat_angle_file(dir, fname);

% Copy mstm files to mstm_fortran for simulation
% source_file = strcat(dir,'/',fname,'.inp');
% source_file = strcat(dir, '/', fname,'_sphere_file.pos');
% source_file = strcat(dir, '/', fname,'_scat_angles.dat');            

% mstm_exe_file = '/home/parkerwray/Desktop/Simulate_Nanoparticle_Clusters/mstm_fortran';
% copyfile(mstm_exe_file, dir);
end
























