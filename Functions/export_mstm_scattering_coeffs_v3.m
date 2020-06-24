
function [sphere, alpha, beta, cbeam, m_med, Nspheres, Nequs] = export_mstm_scattering_coeffs_v3(File)
% clearvars -except File FileName PathName
% clc
%[File, FileName, PathName] = get_file();
line = 1;
fid = fopen(File,'r');
while ~feof(fid)
    text{line} = fgetl(fid);
    line = line+1;
end
%a = fscanf(fid, '%f');
fclose(fid);
text = text.';
%%
sphere = struct();

a = str2num(text{1});
Nspheres = a(1);
Nequs = a(2);
alpha = a(3);
beta = a(4);
cbeam = a(5);
m_med = a(6);



a = str2num(text{2});
sphere(1).Nsphere_order = a(1);
sphere(1).numfieldexp = a(2);
sphere(1).x = a(3);
sphere(1).loc_x = a(4);
sphere(1).loc_y = a(5);
sphere(1).loc_z = a(6);
sphere(1).m_sph_r = a(7);
sphere(1).m_sph_i = a(8);

i = 1;
count = 0;
dummy_a_te = [];
dummy_a_tm = [];
dummy_f_te = [];
dummy_f_tm = [];

for idx = 3:length(text)
    
    a = str2num(text{idx});
    
    if length(a) == 9
        sphere(i).Qext_par = a(1);
        sphere(i).Qext_per = a(2);
        sphere(i).Qabs_par = a(4);
        sphere(i).Qabs_per = a(5);
        sphere(i).Qsca_par = a(7);
        sphere(i).Qsca_per = a(8);
        count = 0;
    end
    if length(a) == 4
        if rem(count,2) == 0
        dummy_a_te = [dummy_a_te, a(1)+1i.*a(2)];
        dummy_f_te = [dummy_f_te, a(3)+1i.*a(4)];  % Actually polarizatin? (perp and para)?
        else 
        dummy_a_tm = [dummy_a_tm, a(1)+1i.*a(2)];
        dummy_f_tm = [dummy_f_tm, a(3)+1i.*a(4)];
        end
        count = count+1;
    end 
    
    if length(a) ==  10
        sphere(i).a_tm_LR = get_coeffs(sphere(i).Nsphere_order, dummy_a_tm);
        sphere(i).a_te_LR = get_coeffs(sphere(i).Nsphere_order, dummy_a_te);
        sphere(i).f_tm_LR = get_coeffs(sphere(i).Nsphere_order, dummy_f_tm);
        sphere(i).f_te_LR = get_coeffs(sphere(i).Nsphere_order, dummy_f_te);
        
        dummy_a_te = [];
        dummy_a_tm = [];
        dummy_f_te = [];
        dummy_f_tm = [];
       
        i = i+1;
        sphere(i).Nsphere_order = a(1);
        sphere(i).numfieldexp = a(2);
        sphere(i).x = a(3);
        sphere(i).loc_x = a(4);
        sphere(i).loc_y = a(5);
        sphere(i).loc_z = a(6);
        sphere(i).m_sph_r = a(7);
        sphere(i).m_sph_i = a(8);
        count = 0;
    end
        
end
sphere(i).a_tm_LR = get_coeffs(sphere(i).Nsphere_order, dummy_a_tm);
sphere(i).a_te_LR = get_coeffs(sphere(i).Nsphere_order, dummy_a_te);

sphere(i).f_tm_LR = get_coeffs(sphere(i).Nsphere_order, dummy_f_tm);
sphere(i).f_te_LR = get_coeffs(sphere(i).Nsphere_order, dummy_f_te);



%% 
%{
The MSTM code works in an LR basis system. 
    L = TE+TM = TE';
    R = TE-TM = TM';
The scattering coefficients are in this basis and need to be converted to
traditional "TE" and "TM" w.r.t. the incident k-vector.

    TE = 0.5.*(TE' + TM');
    TM = 0.5.*(TE' - TM');
%}

% Convert coeffs to incident k-vector TE-TM basis.
for idx = 1:(i)
   for idx2 = 1:sphere(i).Nsphere_order

        sphere(idx).a_te{idx2} = 0.5.*(sphere(idx).a_te_LR{idx2} + sphere(idx).a_tm_LR{idx2}); 
        sphere(idx).a_tm{idx2} = 0.5.*(sphere(idx).a_te_LR{idx2} - sphere(idx).a_tm_LR{idx2}); 
   
   end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTION TO COLLECT COEFFICIENTS BY ORDER. 

%{
A = scattering coefficients
F = plane wave (incident) coefficients

A_mnp = sum_n=1_Norder{  sum_m=-n_n  { sum_p=1+2 {   }}}

n = order
m = degree  (this varies in size based on n)
p = mode (p=1 -> TM, p=2 -> TE)

The cell value is the order (1-Nsphere_order).
Inside each cell are the degrees going from (-n to n)
Modes are seperated by different variables (e.g., a_tm and a_te)

%}

function coeffs = get_coeffs(Norder, coeff_array)

for idx_dummy = 1:Norder
    for i2 = 1:2*idx_dummy+1
        dummy(i2) = coeff_array(i2);
    end
    coeffs{idx_dummy} = dummy;
    coeff_array(1:i2) = [];
end
clearvars idx_dummy

end


end


