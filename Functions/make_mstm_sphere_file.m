
function make_mstm_sphere_file(directory, fname, sph)

inputfilename = strcat(directory, '/', fname,'_sphere_file.pos');
fid = fopen(inputfilename,'w');
for idx = 1:size(sph,1)
    fprintf(fid, [num2str(sph(idx,1)),',',num2str(sph(idx,2)),',',...
    num2str(sph(idx,3)),',',num2str(sph(idx,4)),',',num2str(sph(idx,5))...
    ,',',num2str(sph(idx,6)),'\n']);
end

fclose(fid);
end
   













































