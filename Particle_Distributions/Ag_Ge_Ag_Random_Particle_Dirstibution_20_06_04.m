function spheres = Ag_Ge_Ag_Random_Particle_Dirstibution_20_06_04(wavelengths)

mat_file = '/global/cfs/cdirs/m2542/parkerwray/20_06_15_Material Data/Refractive Index Info/Ag_Yang.csv';
name = 'Ag';
mat1 = pw_read_refinfo_mat(name, mat_file, 1);

mat_file = '/global/cfs/cdirs/m2542/parkerwray/20_06_15_Material Data/Refractive Index Info/Ge.csv';
name = 'Ge';
mat2 = pw_read_refinfo_mat(name, mat_file, 1);
mat2.k(mat2.k==0) = 1.*10^(-6);

%wavelengths = [400:50:500]; %[400:10:1000, 1100:100:5000, 5100:200:10000];
mat1 = index_in_range(mat1, wavelengths, 1);
mat2 = index_in_range(mat2, wavelengths, 1);

%% Compile sphere [radius, x, y, z, n, k] into matrix for outermost shell
clearvars spheres spheres_layer1 spheres_layer2 spheres_layer3

r = 82;
dummy = [r.*ones(size(cord,1),1), cord];
for idx = 1:length(wavelengths)
  spheres_layer3(:,:,idx) = [dummy,...
        mat1.n(idx).*ones(size(cord,1),1),...
        mat1.k(idx).*ones(size(cord,1),1)];
end

r = 80;
dummy = [r.*ones(size(cord,1),1), cord];
for idx = 1:length(wavelengths)
  spheres_layer2(:,:,idx) = [dummy,...
        mat2.n(idx).*ones(size(cord,1),1),...
        mat2.k(idx).*ones(size(cord,1),1)];
end

r = 40;
dummy = [r.*ones(size(cord,1),1), cord];
for idx = 1:length(wavelengths)
  spheres_layer1(:,:,idx) = [dummy,...
        mat1.n(idx).*ones(size(cord,1),1),...
        mat1.k(idx).*ones(size(cord,1),1)];
end

spheres = [spheres_layer3; spheres_layer2; spheres_layer1];


end



















