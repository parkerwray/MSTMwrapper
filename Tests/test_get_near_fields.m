

%% Get near field data
% clear
% clc
% if ~exist('File_NF')
%     [File_NF, FileName, PathName] = get_file();
% end
fid = fopen(File_NF,'r');

line = 1;
while ~feof(fid)
    text{line} = fgetl(fid);
    line = line+1;
end
fclose(fid)
text = text';

line = 1;
dummy = str2num(text{line});
NA = dummy(1);
NB = dummy(2);

line = 2;
Nspheres = str2num(text{line});


for idx = 1:Nspheres
    line = line+1;
    dummy = str2num(text{line});
    
    sphere_kA(idx) = dummy(1);
    sphere_kB(idx) = dummy(2);
    sphere_r(idx) = dummy(3);
end



for idx = 1:length(text)-line
    line = line+1;
    dummy = str2num(text{line});
    kA(idx) = dummy(1);
    kB(idx) = dummy(2);


    Ex(idx) = dummy(3)+1i.*dummy(4);
    Ey(idx) = dummy(5)+1i.*dummy(6);
    Ez(idx) = dummy(7)+1i.*dummy(8);
    Hx(idx) = dummy(9)+1i.*dummy(10);
    Hy(idx) = dummy(11)+1i.*dummy(12);
    Hz(idx) = dummy(13)+1i.*dummy(14);
    
end


%
clearvars E E2
KA = reshape(kA, NB, NA);
KB = reshape(kB, NB, NA);

E(:,:,1) = reshape(Ex, NB, NA);
E(:,:,2) = reshape(Ey, NB, NA);
E(:,:,3) = reshape(Ez, NB, NA);

%
E2 = vecnorm(E,2,3); %Take p=2 norm in 3rd dimension
lda = wavelengths;
k = 2.*pi./lda;
figure, 
hold on 
imagesc((kA./k).*10^(-3), (kB/k).*10^(-3), (E2).^2)
%contourf((KA./k).*10^(-3), (KB./k).*10^(-3), log10(E2), 200, 'edgecolor','none')
circle(0,0,22)
xlabel('Z (um)')
ylabel('X (um)')
colormap('jet')
h = colorbar;
%ylabel(h, 'Log_{10}(|E_T|^2/|E_o|^2)')
ylabel(h, '(|E_T|/|E_o|)')
set(gca, 'FontSize', 24)
pbaspect([1 1 1])










vec = r3-r1-0.2;


function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'w','LineWidth', 2);
end












% %%
% figure, 
% imagesc((kA./k).*10^(-3), (kB/k).*10^(-3), abs(E(:,:,1)))
% %contourf((KA./k).*10^(-3), (KB./k).*10^(-3), log10(E2), 200, 'edgecolor','none')
% xlabel('Y (um)')
% ylabel('Z (um)')
% colormap('jet')
% h = colorbar;
% ylabel(h, '(|E_x|)')
% set(gca, 'FontSize', 24)
% pbaspect([1 1 1])
% 
% %%
% figure, 
% imagesc((kA./k).*10^(-3), (kB/k).*10^(-3), abs(E(:,:,2)))
% %contourf((KA./k).*10^(-3), (KB./k).*10^(-3), log10(E2), 200, 'edgecolor','none')
% xlabel('Y (um)')
% ylabel('Z (um)')
% colormap('jet')
% h = colorbar;
% ylabel(h, '(|E_y|)')
% set(gca, 'FontSize', 24)
% pbaspect([1 1 1])
% %%
% figure, 
% imagesc((kA./k).*10^(-3), (kB/k).*10^(-3), abs(E(:,:,3)))
% %contourf((KA./k).*10^(-3), (KB./k).*10^(-3), log10(E2), 200, 'edgecolor','none')
% xlabel('Y (um)')
% ylabel('Z (um)')
% colormap('jet')
% h = colorbar;
% ylabel(h, '(|E_y|)')
% set(gca, 'FontSize', 24)
% pbaspect([1 1 1])

