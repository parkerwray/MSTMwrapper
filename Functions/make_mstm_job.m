

function make_mstm_job(dir,fname,spheres,mstm_flags,convergence,input_beam,near_field,mstm_input_params) 
            
% Make mstm files to be used
make_mstm_input_file_v3(dir,...
                        fname,...
                        mstm_flags,...
                        convergence,...
                        input_beam,...
                        near_field,...
                        mstm_input_params);
make_mstm_sphere_file(dir, fname, spheres);
make_mstm_scat_angle_file(dir, fname);


end
























