
%{
This code compares mode coefficients of observation particles in random 
particle films to the isolated particle case of same radius.  

Assumed variables existing in the workspace

Modes for isolated particles: Ai Axi Bi Bxi
Modes for comparison particles: A Ax B Bx

FBR for isolated particles: Ipari Iperi FBRi
FBR for comparison particles: Ipar Iper FBR

%} 

clc
r = 140;
ff = 0.2;

[~, idxr] = find(center_radiis == r);
[~, idxri] = find(center_radiis_i == r);
[~, idxff] = find(fill_fractions == ff);
idxn = 1;
idxm = 2;
x2 = (((2.*pi.*r)./wavelengths).^2).';
xi2 = (((2.*pi.*r)./wavelengths_i).^2);

I = squeeze((mean(Ipar,3)+mean(Iper,3))./2);
fbr = (I(:, :, :, 1)./I(:, :, :, end)); 
fbri = (Iipar(:,:,1)./Iipar(:,:,2));

figure, 
yyaxis left
hold on 
% Plot A11
[data_mean, data_std] = get_mag_stat(Astats, idxr, idxff, idxn, idxm);
shadedErrorBar(wavelengths, data_mean./x2, data_std./x2,...
    'lineprops',{'-g','MarkerFaceColor','g', 'LineWidth',6},...
    'patchSaturation', 0.1)
% Plot B11
[data_mean, data_std] = get_mag_stat(Bstats, idxr, idxff, idxn, idxm);
shadedErrorBar(wavelengths, data_mean./x2, data_std./x2,...
    'lineprops',{'-m','MarkerFaceColor','m', 'LineWidth',6},...
    'patchSaturation', 0.1)
%Plot Ax11
[data_mean, data_std] = get_mag_stat(Axstats, idxr, idxff, idxn, idxm);
shadedErrorBar(wavelengths, data_mean./x2, data_std./x2,...
    'lineprops',{'-c','MarkerFaceColor','c', 'LineWidth',6},...
    'patchSaturation', 0.1)
% Plot Bx11
[data_mean, data_std] = get_mag_stat(Bxstats, idxr, idxff, idxn, idxm);
shadedErrorBar(wavelengths, data_mean./x2, data_std./x2,...
    'lineprops',{'-r','MarkerFaceColor','r', 'LineWidth',6},...
    'patchSaturation', 0.1)
% Plot A21
[data_mean, data_std] = get_mag_stat(Astats, idxr, idxff, idxn+1, idxm);
shadedErrorBar(wavelengths, data_mean./x2, data_std./x2,...
    'lineprops',{'-b','MarkerFaceColor','b', 'LineWidth',6},...
    'patchSaturation', 0.1)
% Plot B21
[data_mean, data_std] = get_mag_stat(Bstats, idxr, idxff, idxn+1, idxm);
shadedErrorBar(wavelengths, data_mean./x2, data_std./x2,...
    'lineprops',{'-y','MarkerFaceColor','y', 'LineWidth',6},...
    'patchSaturation', 0.1)

% Plot isolated A11
plot(wavelengths_i, (abs(squeeze(Ai(idxri,:, idxn, idxm))).^2)./xi2, ...
    ':g', 'LineWidth', 8)
% Plot isolated B11
plot(wavelengths_i, (abs(squeeze(Bi(idxri,:, idxn, idxm))).^2)./xi2, ...
    ':m', 'LineWidth', 8)
% Plot isolated A21
plot(wavelengths_i, (abs(squeeze(Ai(idxri,:, idxn+1, idxm))).^2)./xi2, ...
    ':b', 'LineWidth', 8)
% Plot isolated B21
plot(wavelengths_i, (abs(squeeze(Bi(idxri,:, idxn+1, idxm))).^2)./xi2, ...
    ':y', 'LineWidth', 8)

ylim([0,1])
ylabel('Efficiency')

yyaxis right
plot(wavelengths, log10(squeeze(fbr(idxr, idxff, :))),'k', 'LineWidth', 6)
plot(wavelengths_i, log10(squeeze(fbri(idxri,:))), ':k', 'LineWidth', 8)
ylabel('FBR')

hold off
title(['Random distribution versus isolated particle', newline,...
    'Normal distribution radius =  ', num2str(r_mean),'+/-',num2str(r_sigma),'nm',newline,...
    'Observation particle radius = ', num2str(r), 'nm with fill fraction = ', num2str(ff*100),'%'])
legend('A_{11}','B_{11}','Ax_{11}','Bx_{11}','A_{21}','B_{21}',...
    'a_{1}','b_{1}','a_{2}','b_{2}'),
xlim([500,900])
xlabel('Wavelength (nm)')



%%

figure, 
hold on 
% Plot A11
[data_mean, data_std] = get_phase_stat(Astats, idxr, idxff, idxn, idxm);
shadedErrorBar(wavelengths, data_mean, data_std,...
    'lineprops',{'-g','MarkerFaceColor','g', 'LineWidth',6},...
    'patchSaturation', 0.1)
% Plot B11
[data_mean, data_std] = get_phase_stat(Bstats, idxr, idxff, idxn, idxm);
shadedErrorBar(wavelengths, data_mean, data_std,...
    'lineprops',{'-m','MarkerFaceColor','m', 'LineWidth',6},...
    'patchSaturation', 0.1)
%Plot Ax11
[data_mean, data_std] = get_phase_stat(Axstats, idxr, idxff, idxn, idxm);
shadedErrorBar(wavelengths, data_mean, data_std,...
    'lineprops',{'-c','MarkerFaceColor','c', 'LineWidth',6},...
    'patchSaturation', 0.1)
% Plot Bx11
[data_mean, data_std] = get_phase_stat(Bxstats, idxr, idxff, idxn, idxm);
shadedErrorBar(wavelengths, data_mean, data_std,...
    'lineprops',{'-r','MarkerFaceColor','r', 'LineWidth',6},...
    'patchSaturation', 0.1)
% Plot A21
[data_mean, data_std] = get_phase_stat(Astats, idxr, idxff, idxn+1, idxm);
shadedErrorBar(wavelengths, data_mean, data_std,...
    'lineprops',{'-b','MarkerFaceColor','b', 'LineWidth',6},...
    'patchSaturation', 0.1)
% Plot B21
[data_mean, data_std] = get_phase_stat(Bstats, idxr, idxff, idxn+1, idxm);
shadedErrorBar(wavelengths, data_mean, data_std,...
    'lineprops',{'-y','MarkerFaceColor','y', 'LineWidth',6},...
    'patchSaturation', 0.1)

% Plot isolated A11
plot(wavelengths_i, (angle(squeeze(Ai(idxri,:, idxn, idxm))).*180/pi    ), ...
    ':g', 'LineWidth', 8)
% Plot isolated B11
plot(wavelengths_i, (angle(squeeze(Bi(idxri,:, idxn, idxm))).*180/pi  ), ...
    ':m', 'LineWidth', 8)
% Plot isolated A21
plot(wavelengths_i, (angle(squeeze(Ai(idxri,:, idxn+1, idxm))).*180/pi  ), ...
    ':b', 'LineWidth', 8)
% Plot isolated B21
plot(wavelengths_i, (angle(squeeze(Bi(idxri,:, idxn+1, idxm))).*180/pi  ), ...
    ':y', 'LineWidth', 8)

hold off
title(['Random distribution versus isolated particle', newline,...
    'Normal distribution radius =  ', num2str(r_mean),'+/-',num2str(r_sigma),'nm',newline,...
    'Observation particle radius = ', num2str(r), 'nm with fill fraction = ', num2str(ff*100),'%'])
legend('A_{11}','B_{11}','Ax_{11}','Bx_{11}','A_{21}','B_{21}',...
    'a_{1}','b_{1}','a_{2}','b_{2}'),
xlim([500,900])
% ylim([0,1])
xlabel('Wavelength (nm)')
ylabel('Phase')





























function [data_mean, data_std] = get_mag_stat(data, idxr, idxff, idxn, idxm)

    data_mean = squeeze(data.mean.mag(idxr,idxff,1,:,idxn,idxm));
    data_std = squeeze(data.std.mag(idxr,idxff,1,:,idxn,idxm));
    
end


function [data_mean, data_std] = get_phase_stat(data, idxr, idxff, idxn, idxm)

    data_mean = squeeze(data.mean.phase(idxr,idxff,1,:,idxn,idxm));
    data_std = squeeze(data.std.phase(idxr,idxff,1,:,idxn,idxm));
    
end




























