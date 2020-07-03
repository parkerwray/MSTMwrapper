function spheres = give_material_properties(spheres, mat_file, wavelength)




% Give spheres material properties
%mat_file = '/global/cfs/cdirs/m2542/parkerwray/20_06_15_Material Data/Refractive Index Info/Ag_Yang.csv';
name = 'name';
mat = index_in_range(pw_read_refinfo_mat(name, mat_file, 0), wavelength, 0);


spheres = [spheres,...
    mat.n.*ones(size(spheres,1),1),...
    mat.k.*ones(size(spheres,1),1)];








end





