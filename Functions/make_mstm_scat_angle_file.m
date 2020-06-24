
function make_mstm_scat_angle_file(directory, fname)

inputfilename = strcat(directory, '/', fname,'_scat_angles.dat');
fid = fopen(inputfilename,'w');
ang = 0:1:359;
for idx = 1:length(ang)
    fprintf(fid, [num2str(ang(idx)),',',num2str(0),'\n']);
end

fclose(fid);
end