
function [cluster, sphere, excitation, FLAG] = export_output_file_v3(File)

%{
This function reads the standard output file from MSTM and saves the results
into two structure. 

1) sphere structure [array]: Each individual sphere is given a structure 
that houses its information about its size parameter, index, spatial 
coordinates as well as cross sections (extinction, absorption, and 
scattering). The structure also houses the spheres coefficients. These 
coefficients are saved in another text file and grabbed by another 
sub-function.  

2) Cluster structure: This structure houses information related to the
cluser as a whole. This includes the "equivalent sphere" size parameter,
total cluster cross sections (extinction, absorption, adn scattering), 
total cluster asymetery parameter, and the total cluster Mueller Matrix. 

3) 


%}

FLAG = 0;
fid = fopen(File,'r');
if fid == -1
    disp('Output file does not exist!')
    %keyboard;
    cluster = [NaN];
    sphere = [NaN];
    excitation = [NaN];
    keyboard;
    FLAG = 1;
    return;
end

sphere = struct();
cluster = struct();
excitation = struct();
S = struct();

line = 1;
while ~feof(fid)
    text{line} = fgetl(fid);
    line = line+1;
end
fclose(fid);
text = text';


description = ' incident azimuth and polar angles:';
line_offset = 1;
data = read_number_line(text, description, line_offset);
excitation.azimuth = data(:,1);
excitation.polar = data(:,2);

description = ' scattering plane angle:';
line_offset = 1;
data = read_number_line(text, description, line_offset);
excitation.scat_plane_angle = data(:,1);


%
description =  'number of spheres, volume size parameter:';
line_offset = 1;
data = read_number_line(text, description, line_offset);

cluster.Nspheres = data(:,1);
cluster.ka = data(:,2);

%
for i = 1:cluster.Nspheres(1)
    description =  ' sphere host   ka     x-x(host) y-y(host) z-z(host)  Re(m)     Im(m)      Qext         Qsca        Qabs         Qabs(V)';
    data = read_number_line(text, description, i);

    sphere(i).ka = data(:,3);
    sphere(i).cords  = data(:,4:6);
    ri = data(:,7:8);
    sphere(i).m = ri(1)+1i.*ri(2);
    sphere(i).qe = data(:,9);
    sphere(i).qs = data(:,10);
    sphere(i).qa = data(:,11);
    sphere(i).qav = data(:,12);

end

%
description =  'total ext, abs, scat efficiencies, w.r.t. xv, and asym. parm';
line_offset = 1;
data = read_number_line(text, description, line_offset);
cluster.qext = data(:,1);
cluster.qabs = data(:,2);
cluster.qsca = data(:,3);
cluster.asym = data(:,4);

%%
% Get elements of the phase matrix (the matrix that describes the stokes
% parameters) for the entire cluster as an effectie particle. 
description =  'number scattering angles:';
line_offset = 1;
data = read_number_line(text, description, line_offset);
cluster.Ntheta = data(:,1);
description =  ' scattering matrix elements';
elements = read_text_line(text, description, 1);

for i = 1:cluster.Ntheta(1)

    data = read_number_line(text, elements, i);
    
    c = 1;
    S(i).theta = data(:,c);
    
    if strfind(elements,'phi')
        c = c+1;
        S(i).phi = data(:,c);
    else
        S(i).phi = 0;        
    end
    
    if strfind(elements,'11')
        c = c+1;
        S(i).s11 = data(:,c);
    else
        S(i).s11 = 0;        
    end
    
    if strfind(elements,'12') 
        c = c+1;
        S(i).s12 = data(:,c);  
     else
        S(i).s12 = 0;         
    end
    
    if strfind(elements,'13')    
        c = c+1;
        S(i).s13 = data(:,c);
    else
        S(i).s13 = 0;          
    end
    
    if strfind(elements,'14')    
        c = c+1;
        S(i).s14 = data(:,c); 
    else
        S(i).s14 = 0;          
    end
    
    if strfind(elements,'21')    
        c = c+1;
        S(i).s21 = data(:,c); 
    else
        S(i).s21 = 0;          
    end
    
    if strfind(elements,'22')  
        c = c+1;
        S(i).s22 = data(:,c); 
    else
        S(i).s22 = 0;          
    end
    
    if strfind(elements,'23')    
        c = c+1;
        S(i).s23 = data(:,c); 
    else
        S(i).s23 = 0;          
    end    
    
    if strfind(elements,'24')  
        c = c+1;
        S(i).s24 = data(:,c);  
    else
        S(i).s24 = 0;         
    end    
    
    if strfind(elements,'31')   
        c = c+1;
        S(i).s31 = data(:,c); 
    else
        S(i).s31 = 0;          
    end
    
    if strfind(elements,'32')  
        c = c+1;
        S(i).s32 = data(:,c);  
    else
        S(i).s32 = 0;          
    end
    
    if strfind(elements,'33')   
        c = c+1;
        S(i).s33 = data(:,c); 
     else
        S(i).s33 = 0;         
    end
    
    if strfind(elements,'34')       
        c = c+1;
        S(i).s34 = data(:,c);  
    else
        S(i).s34 = 0;          
    end
    
    if strfind(elements,'41')   
        c = c+1;
        S(i).s41 = data(:,c); 
    else
        S(i).s41 = 0;          
    end
    
    if strfind(elements,'42')  
        c = c+1;
        S(i).s42 = data(:,c);  
    else
        S(i).s42 = 0;          
    end
    
    if strfind(elements,'43')   
        c = c+1;
        S(i).s43 = data(:,c);
    else
        S(i).s43 = 0;          
    end
    
    if strfind(elements,'44')       
        c = c+1;
        S(i).s44 = data(:,c);  
    else
        S(i).s44 = 0;          
    end
    
end

cluster.S = S;
clearvars S


%% From the phase matrix (Z) derive the amplitude matrix (S) 

% % % rho = (1/sqrt(2)).*[1,0,0,1; 1,0,0,-1; 0,-1,-1,0; 0,-i,i,0];
% % % rho_inv = (1/sqrt(2)).*[1,1,0,0; 0,0,-1,i; 0,0,-1,-i; 1,-1,0,0];
% % % 
% % % 
% % % 
% % % S_kron = rho_inv*cluster.Z*rho;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function line_below = read_number_line(text, description, line_offset)

    
    idx = strfind(text,description);
    counter = 1;
    for idx_dummy = 1:length(text)
        if idx{idx_dummy}>0
            line_below{counter} = str2num(text{idx_dummy+line_offset});
            counter = counter+1;
        end
    end
    line_below = cell2mat(line_below.');
    clearvars idx idx_dummy
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function line_below = read_text_line(text, description, line_offset)

    
    idx = strfind(text,description);

    counter = 1;
    for idx_dummy = 1:length(text)
        if idx{idx_dummy}>0
            line_below{counter} = text{idx_dummy+line_offset};
            counter = counter+1;
        end
    end
    line_below = cell2mat(line_below.');
    clearvars idx idx_dummy   
end








end

