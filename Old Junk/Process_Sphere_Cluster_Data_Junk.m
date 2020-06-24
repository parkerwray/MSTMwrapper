%% 
%{
This code makes figures for the Huygens monolayer nanoparticle films. 
It calculates the response of a nanoparticle to random particles around it,
then makes relevant figures to show Huygens is preserved, on average. 
%}

%% 

clc
[modes_perp, modes_par, modes] = grab_modes(sphere_coeffs_perp, sphere_coeffs_par);
[A11_mag, B11_mag, Ax11_mag, Bx11_mag] = get_AB_mag_stats(modes);
[I0_par_c, I0_perp_c, I180_par_c, I180_perp_c] = get_I_compact(modes);
[qe,qs,qa, qds] = get_efficiencies(sphere_result_par);





%%
clc
















%%
figure, 
plot(ff_calc)

title('calculated fill fraction')



%%

fig = figure; 
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
s = shadedErrorBar(wavelengths, 6.*A11_mag(:,1)./(modes.x(:,1).^2), 6.*A11_mag(:,2)./(modes.x(:,1).^2),...
    'lineprops',{'-g','MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)


s = shadedErrorBar(wavelengths, 6.*B11_mag(:,1)./(modes.x(:,1).^2), 6.*B11_mag(:,2)./(modes.x(:,1).^2),...
    'lineprops',{'-m','MarkerFaceColor','m', 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)


s = shadedErrorBar(wavelengths, 6.*Ax11_mag(:,1)./(modes.x(:,1).^2), 6.*Ax11_mag(:,2)./(modes.x(:,1).^2),...
    'lineprops',{'-r','MarkerFaceColor','r', 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)


s = shadedErrorBar(wavelengths, 6.*Bx11_mag(:,1)./(modes.x(:,1).^2), 6.*Bx11_mag(:,2)./(modes.x(:,1).^2),...
    'lineprops',{'-c','MarkerFaceColor','c', 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)


%ylabel('Efficiency ( $$ \frac{6}{x^{2}}E[||\dot||^{2}] $$) ','interpreter','latex')
ylabel('Efficiency (E[\bullet]/x^{2}) ' )
%ylim([0,2.2])
yyaxis right
plot(wavelengths, log10(mean(I0_perp_c,2)./mean(I180_perp_c,2)), 'LineWidth',6)
ylabel('Log_{10} forward-to-backward ratio')
ylim([min(log10(mean(I0_perp_c,2)./mean(I180_perp_c,2)))-0.1, max(log10(mean(I0_perp_c,2)./mean(I180_perp_c,2)))+0.1])
%ylim([-0.5,1.5])
legend('$$ ||A_{11}||^{2}$$','$$  ||B_{11}||^{2}$$',...
    '$$   ||A_{11}^{\times}||^{2}$$','$$ ||B_{11}^{\times}||^{2}$$',...
    'FBR',...
    'interpreter','latex');

%xlim([250,600])
%xticks(375:25:500)
title('30% Fill fraction','FontSize', 34)
xlabel('Wavelength (nm)')

pbaspect([1 1 1])
set(gca,'FontSize',34)

%%
figure, plot(wavelengths, mean(qs+qds,2))

%%
clc
q = mean((qs+qds),2);
%q = interp1(wavelengths, q , lda_sim);
% max_q = 100/(Nspheres*pi);
ffq = ff.*q;
% ffq(ffq>max_q) = max_q;

fbr = mean(I0_perp_c,2)./mean(I180_perp_c,2);
%fbr = interp1(wavelengths, fbr, lda_sim);
% Rlum = interp1(lda, R20m, wavelengths).';
% Alum = interp1(lda, A20m, wavelengths).';

%A = Alum;
A = ff.*mean(qa,2);

R = ffq./(1+fbr);
figure, 

plot(wavelengths, 100.*R, 'b', 'LineWidth', 6)
hold on 
plot(wavelengths, 100.*(A), 'r', 'LineWidth', 6)
plot(wavelengths, 100.*(1-A-R), 'k', 'LineWidth', 6)
%plot(lda, 100.*T20m, 'k', 'LineWidth', 6)
hold off

pbaspect([1 1 1])
set(gca,'FontSize',34)
legend('Reflection (calc.)', ...
    'Absorption (calc.)', 'Transmission (calc.)','FontSize', 34)


%xlim([500,850])
title('10% Fill fraction','FontSize', 34)
xlabel('Wavelength (nm)','FontSize', 34)
ylabel('Avg. transmission or reflection (%)','FontSize', 34)




%ylim([0,0.4])
% yyaxis right
% plot(lda_sim, factor)





%%
clc
q = mean((qs+qds),2);
q = interp1(wavelengths, q , lda_sim);
max_q = 400/(19*4*pi);
q( q > max_q) = max_q;


q = mean((qe),2);
q = interp1(wavelengths, q , lda_sim);




fbr = mean(I0_perp_c,2)./mean(I180_perp_c,2);
fbr = interp1(wavelengths, fbr, lda_sim);
I0 = interp1(wavelengths, mean(I0_perp_c,2)./((x.').^2), lda_sim);

T = 1+(I0-q).*0.6;

figure, 
plot(lda_sim, mat2gray(T))
hold on 
plot(lda_sim, mat2gray(T60m))
hold off

% R = (1-A60m./(q./max_q))./(1+fbr);
% 
% 
% figure, 
% 
% plot(lda_sim, R)
% hold on 
% % shadedErrorBar(lda_sim, R60m, R60s,...
% %     'lineprops',{'-c','MarkerFaceColor','c', 'LineWidth',6},'patchSaturation', 0.1)
% plot(lda_sim, R60m)
% hold off
% legend('test', 'actual')
% ylim([0,0.4])
% % yyaxis right
% % plot(lda_sim, factor)
% 
% 
% figure, 
% plot(lda_sim, (1-A60m-R))
% hold on 
% plot(lda_sim, T60m)
% hold off


% 
% figure, plot(wavelengths, mean(IA0_perp_c,2))

%% Plot Higher Order Mode Magnitudes
clc
%wavelengths = 300:2:500;
x = mean(modes.x,2).';
Amn_mag_m = squeeze(mean(abs(modes.Amn).^2,3));
Amn_mag_s = squeeze(std(abs(modes.Amn).^2,[],3));
Bmn_mag_m = squeeze(mean(abs(modes.Bmn).^2,3));
Bmn_mag_s = squeeze(std(abs(modes.Bmn).^2,[],3));

Axmn_mag_m = squeeze(mean(abs(modes.Axmn).^2,3));
Axmn_mag_s = squeeze(std(abs(modes.Axmn).^2,[],3));
Bxmn_mag_m = squeeze(mean(abs(modes.Bxmn).^2,3));
Bxmn_mag_s = squeeze(std(abs(modes.Bxmn).^2,[],3));




figure, 
s = shadedErrorBar(wavelengths, 6.*Amn_mag_m(1,:)./(x.^2), 0.*Amn_mag_s(1,:)./(x.^2),...
    'lineprops',{'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
s = shadedErrorBar(wavelengths, 10.*Amn_mag_m(3,:)./(x.^2), 0.*Amn_mag_s(3,:)./(x.^2),...
    'lineprops',{'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
s = shadedErrorBar(wavelengths, 10.*Amn_mag_m(4,:)./(x.^2), 0.*Amn_mag_s(4,:)./(x.^2),...
    'lineprops',{'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
s = shadedErrorBar(wavelengths, 10.*Amn_mag_m(5,:)./(x.^2), 0.*Amn_mag_s(5,:)./(x.^2),...
    'lineprops',{'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
s = shadedErrorBar(wavelengths, 6.*Bmn_mag_m(1,:)./(x.^2), 0.*Bmn_mag_s(1,:)./(x.^2),...
    'lineprops',{ 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
s = shadedErrorBar(wavelengths, 10.*Bmn_mag_m(3,:)./(x.^2), 0.*Bmn_mag_s(3,:)./(x.^2),...
    'lineprops',{'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
s = shadedErrorBar(wavelengths, 10.*Bmn_mag_m(4,:)./(x.^2), 0.*Bmn_mag_s(4,:)./(x.^2),...
    'lineprops',{'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
s = shadedErrorBar(wavelengths, 10.*Bmn_mag_m(5,:)./(x.^2), 0.*Bmn_mag_s(5,:)./(x.^2),...
    'lineprops',{'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
% shadedErrorBar(wavelengths, 6.*Axmn_mag_m(1,:)./(x.^2), 0.*Axmn_mag_s(1,:)./(x.^2),...
%     'lineprops',{':', 'LineWidth',6},'patchSaturation', 0.1)

% shadedErrorBar(wavelengths, 10.*Axmn_mag_m(3,:)./(x.^2), 0.*Axmn_mag_s(3,:)./(x.^2),...
%     'lineprops',{':', 'LineWidth',6},'patchSaturation', 0.1)

s = shadedErrorBar(wavelengths, 10.*Axmn_mag_m(4,:)./(x.^2), 0.*Axmn_mag_s(4,:)./(x.^2),...
    'lineprops',{':', 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
s = shadedErrorBar(wavelengths, 10.*Axmn_mag_m(5,:)./(x.^2), 0.*Axmn_mag_s(5,:)./(x.^2),...
    'lineprops',{':', 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
% shadedErrorBar(wavelengths, 6.*Bxmn_mag_m(1,:)./(x.^2), 0.*Bxmn_mag_s(1,:)./(x.^2),...
%     'lineprops',{':',  'LineWidth',6},'patchSaturation', 0.1)

% shadedErrorBar(wavelengths, 10.*Bxmn_mag_m(3,:)./(x.^2), 0.*Bxmn_mag_s(3,:)./(x.^2),...
%     'lineprops',{':', 'LineWidth',6},'patchSaturation', 0.1)

s = shadedErrorBar(wavelengths, 10.*Bxmn_mag_m(4,:)./(x.^2), 0.*Bxmn_mag_s(4,:)./(x.^2),...
    'lineprops',{':', 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
s = shadedErrorBar(wavelengths, 10.*Bxmn_mag_m(5,:)./(x.^2), 0.*Bxmn_mag_s(5,:)./(x.^2),...
    'lineprops',{':', 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)
set(gca, 'YScale', 'log')
% ylim([0,0.06])
ylabel('Efficiency (E[\bullet]/x^{2}) ' )
% legend('$$ ||A_{01}||^{2}$$','$$  ||A_{02}||^{2}$$',...
%     '$$   ||A_{12}||^{2}$$','$$ ||A_{22}||^{2}$$',...
%     '$$ ||B_{01}||^{2}$$','$$  ||B_{02}||^{2}$$',...
%     '$$   ||B_{12}||^{2}$$','$$ ||B_{22}||^{2}$$',...
%     '$$ ||A^{\times}_{01}||^{2}$$','$$  ||A^{\times}_{02}||^{2}$$',...
%     '$$   ||A^{\times}_{12}||^{2}$$','$$ ||A^{\times}_{22}||^{2}$$',...
%     '$$ ||B^{\times}_{01}||^{2}$$','$$  ||B^{\times}_{02}||^{2}$$',...
%     '$$   ||B^{\times}_{12}||^{2}$$','$$ ||B^{\times}_{22}||^{2}$$',...
%     'interpreter','latex');

legend('$$ ||A_{01}||^{2}$$','$$  ||A_{02}||^{2}$$',...
    '$$   ||A_{12}||^{2}$$','$$ ||A_{22}||^{2}$$',...
    '$$ ||B_{01}||^{2}$$','$$  ||B_{02}||^{2}$$',...
    '$$   ||B_{12}||^{2}$$','$$ ||B_{22}||^{2}$$',...
    '$$   ||A^{\times}_{12}||^{2}$$','$$ ||A^{\times}_{22}||^{2}$$',...
    '$$   ||B^{\times}_{12}||^{2}$$','$$ ||B^{\times}_{22}||^{2}$$',...
    'interpreter','latex');

box on 
ylim([10.^(-10),1])
%xlim([375,500])
%xticks(375:25:500)
title('30% Fill fraction','FontSize', 34)
xlabel('Wavelength (nm)')

pbaspect([1 1 1])
set(gca,'FontSize',34)


%% Plot 1st order mode phases

clc
%wavelengths = 350:2:550;
x = mean(modes.x,2).';
A11_phase_m1 = (180./pi).*squeeze(mean(wrapTo2Pi(angle(modes.Amn)),3))-360;
A11_phase_s1 = (180./pi).*squeeze(std(wrapTo2Pi(angle(modes.Amn)),[],3));

A11_phase_m2 = (180./pi).*squeeze(mean(wrapToPi(angle(modes.Amn)),3));
A11_phase_s2 = (180./pi).*squeeze(std(wrapToPi(angle(modes.Amn)),[],3));


B11_phase_m = (180./pi).*squeeze(mean(wrapToPi(angle(modes.Bmn)),3));
B11_phase_s = (180./pi).*squeeze(std(wrapToPi(angle(modes.Bmn)),[],3));

Ax11_phase_m = (180./pi).*squeeze(mean(wrapToPi(angle(modes.Axmn)),3));
Ax11_phase_s = (180./pi).*squeeze(std(wrapToPi(angle(modes.Axmn)),[],3));
Bx11_phase_m = (180./pi).*squeeze(mean(wrapToPi(angle(modes.Bxmn)),3));
Bx11_phase_s = (180./pi).*squeeze(std(wrapToPi(angle(modes.Bxmn)),[],3));

% A11_phase_m = [A11_phase_m2(:,1:50),A11_phase_m1(:,51:91)];
% A11_phase_s = [A11_phase_s2(:,1:50),A11_phase_s1(:,51:91)];
figure, 
plot(A11_phase_m1(2,:))
hold on 
plot(A11_phase_m2(2,:))
%plot(A11_phase_m(2,:))
hold off
legend('Option1','Option2', 'Combined')
%%
fig = figure; 
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
s = shadedErrorBar(wavelengths, A11_phase_m(2,:), A11_phase_s(2,:),...
    'lineprops',{'-g','MarkerFaceColor','g', 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)

s = shadedErrorBar(wavelengths, B11_phase_m(2,:), B11_phase_s(2,:),...
    'lineprops',{'-m','MarkerFaceColor','m', 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)

s = shadedErrorBar(wavelengths, Ax11_phase_m(2,:), Ax11_phase_s(2,:),...
    'lineprops',{'-r','MarkerFaceColor','r', 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)

s = shadedErrorBar(wavelengths, Bx11_phase_m(2,:), Bx11_phase_s(2,:),...
    'lineprops',{'-c','MarkerFaceColor','c', 'LineWidth',6},'patchSaturation', 0.1);
set(s.edge, 'LineWidth',0.5)

%ylabel('Efficiency ( $$ \frac{6}{x^{2}}E[||\dot||^{2}] $$) ','interpreter','latex')
ylabel('Phase (degree) ' )
legend('$$ \angle A_{11}$$','$$  \angle B_{11}$$',...
    '$$   \angle A_{11}^{\times}$$','$$ \angle B_{11}^{\times}$$',...
    'FBR',...
    'interpreter','latex');

box on 
ylim([-200,150])
xlim([375,500])
xticks(375:25:500)
title('30% Fill fraction','FontSize', 34)
xlabel('Wavelength (nm)')

pbaspect([1 1 1])
set(gca,'FontSize',34)




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function plot_1st_order_modes_eff(wavelengths, modes_perp, modes_par, modes)


a111_mag_perp = abs(modes_perp.anm1(2,:,:)).^2; %./(modes.x(:,1).^2);
a112_mag_perp = abs(modes_perp.anm2(2,:,:)).^2; %./(modes.x(:,1).^2);
b111_mag_perp = abs(modes_perp.bnm1(2,:,:)).^2; %./(modes.x(:,1).^2);
b112_mag_perp = abs(modes_perp.bnm2(2,:,:)).^2; %./(modes.x(:,1).^2);

a111_mag_par = abs(modes_par.anm1(2,:,:)).^2; %./(modes.x(:,1).^2);
a112_mag_par = abs(modes_par.anm2(2,:,:)).^2; %./(modes.x(:,1).^2);
b111_mag_par = abs(modes_par.bnm1(2,:,:)).^2; %./(modes.x(:,1).^2);
b112_mag_par = abs(modes_par.bnm2(2,:,:)).^2; %./(modes.x(:,1).^2);

A11_ang = (180./pi).*(squeeze(     mean(squeeze(angle(modes.Amn(2,:,:))),2)      ).^2); %./(modes.x(:,1).^2);
B11_ang = (180./pi).*(squeeze(    mean(squeeze(angle(modes.Bmn(2,:,:))),2)      ).^2); %./(modes.x(:,1).^2);
Ax11_ang = (180./pi).*(squeeze(    mean(squeeze(angle(modes.Axmn(2,:,:))),2)      ).^2); %./(modes.x(:,1).^2);
Bx11_ang = (180./pi).*(squeeze(    mean(squeeze(angle(modes.Bxmn(2,:,:))),2)      ).^2); %./(modes.x(:,1).^2);

a111_ang_perp = (180./pi).*angle(modes_perp.anm1(2,:,:)).^2; %./(modes.x(:,1).^2);
a112_ang_perp = (180./pi).*angle(modes_perp.anm2(2,:,:)).^2; %./(modes.x(:,1).^2);
b111_ang_perp = (180./pi).*angle(modes_perp.bnm1(2,:,:)).^2; %./(modes.x(:,1).^2);
b112_ang_perp = (180./pi).*angle(modes_perp.bnm2(2,:,:)).^2; %./(modes.x(:,1).^2);

a111_ang_par = (180./pi).*angle(modes_par.anm1(2,:,:)).^2; %./(modes.x(:,1).^2);
a112_ang_par = (180./pi).*angle(modes_par.anm2(2,:,:)).^2; %./(modes.x(:,1).^2);
b111_ang_par = (180./pi).*angle(modes_par.bnm1(2,:,:)).^2; %./(modes.x(:,1).^2);
b112_ang_par = (180./pi).*angle(modes_par.bnm2(2,:,:)).^2; %./(modes.x(:,1).^2);


[I0_par, I0_perp, I180_par, I180_perp] = get_I(modes_perp, modes_par);
[I0_par_c, I0_perp_c, I180_par_c, I180_perp_c] = get_I_compact(modes);

figure, 
plot(wavelengths, squeeze(abs(modes.Amn(2,:,1))))
hold on 
plot(wavelengths, squeeze(abs(modes.Amn(2,:,2))))
hold off
legend('Anm first', 'Anm second')

figure, 
plot(wavelengths, squeeze(abs(modes.Bmn(2,:,1))))
hold on 
plot(wavelengths, squeeze(abs(modes.Bmn(2,:,2))))
hold off
legend('Bnm first', 'Bnm second')

figure,
plot(wavelengths, A11_mag, '-','LineWidth',8)
hold on 
plot(wavelengths, a112_mag_perp, '-.','LineWidth',8)
plot(wavelengths, a111_mag_par, '-.','LineWidth',8)
plot(wavelengths, B11_mag, '-','LineWidth',6)
plot(wavelengths, b111_mag_perp, '-.','LineWidth',6)
plot(wavelengths, b112_mag_par, '-.','LineWidth',6)
hold off
title('primary modes')
legend('A11','a112perp','a111par',...
        'B11','b111perp','b112par');

figure,
plot(wavelengths, Ax11_mag, '-','LineWidth',8)
hold on 
plot(wavelengths, a111_mag_perp, '-.','LineWidth',8)
plot(wavelengths, a112_mag_par, '-.','LineWidth',8)
plot(wavelengths, Bx11_mag, '-','LineWidth',6)
plot(wavelengths, b112_mag_perp, '-.','LineWidth',6)
plot(wavelengths, b111_mag_par, '-.','LineWidth',6)
hold off
title('conversion modes')
legend('Ax11','a111perp','a112par',...
        'Bx11','b112perp','b111par');        

figure,
plot(wavelengths, A11_ang, '-','LineWidth',8)
hold on 
plot(wavelengths, a112_ang_perp, '-.','LineWidth',8)
plot(wavelengths, a111_ang_par, '-.','LineWidth',8)
plot(wavelengths, B11_ang, '-','LineWidth',6)
plot(wavelengths, b111_ang_perp, '-.','LineWidth',6)
plot(wavelengths, b112_ang_par, '-.','LineWidth',6)
hold off
title('primary modes')
legend('A11','a112perp','a111par',...
        'B11','b111perp','b112par');

figure,
plot(wavelengths, Ax11_ang, '-','LineWidth',8)
hold on 
plot(wavelengths, a111_ang_perp, '-.','LineWidth',8)
plot(wavelengths, a112_ang_par, '-.','LineWidth',8)
plot(wavelengths, Bx11_ang, '-','LineWidth',6)
plot(wavelengths, b112_ang_perp, '-.','LineWidth',6)
plot(wavelengths, b111_ang_par, '-.','LineWidth',6)
hold off
title('conversion modes')
legend('Ax11','a111perp','a112par',...
        'Bx11','b112perp','b111par');   
    
    
    
figure, 
plot(wavelengths, log10(I0_perp./I180_perp),'-.','LineWidth',6)
hold on 
plot(wavelengths, log10(I0_par./I180_par),'-.','LineWidth',8)  
plot(wavelengths, log10(mean(I0_perp_c,2)./mean(I180_perp_c,2)),'-','LineWidth',6)      
plot(wavelengths, log10(mean(I0_par_c,2)./mean(I180_par_c,2)),':','LineWidth',6) 

plot(wavelengths, log10((I0_perp+I0_par)./(I180_perp+I180_par)),'-.','LineWidth',4)
hold off
legend('FBRperp','FBRpar', 'FBRperp C', 'FBRpar C','FBR perp par avg')

%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A11_mag, B11_mag, Ax11_mag, Bx11_mag] = get_AB_mag_stats(modes)

A11_mag_avg = (squeeze(  mean(squeeze(abs(modes.Amn(2,:,:)).^2),2)  )); %./(modes.x(:,1).^2);
B11_mag_avg = (squeeze(  mean(squeeze(abs(modes.Bmn(2,:,:)).^2),2)  )); %./(modes.x(:,1).^2);
Ax11_mag_avg = (squeeze( mean(squeeze(abs(modes.Axmn(2,:,:)).^2),2) )); %./(modes.x(:,1).^2);
Bx11_mag_avg = (squeeze( mean(squeeze(abs(modes.Bxmn(2,:,:)).^2),2) )); %./(modes.x(:,1).^2);


A11_mag_std = (squeeze(  std(squeeze(abs(modes.Amn(2,:,:)).^2),[],2)  )); %./(modes.x(:,1).^2);
B11_mag_std = (squeeze(  std(squeeze(abs(modes.Bmn(2,:,:)).^2),[],2)  )); %./(modes.x(:,1).^2);
Ax11_mag_std = (squeeze( std(squeeze(abs(modes.Axmn(2,:,:)).^2),[],2) )); %./(modes.x(:,1).^2);
Bx11_mag_std = (squeeze( std(squeeze(abs(modes.Bxmn(2,:,:)).^2),[],2) )); %./(modes.x(:,1).^2);

A11_mag = [A11_mag_avg,A11_mag_std];
B11_mag = [B11_mag_avg,B11_mag_std];
Ax11_mag = [Ax11_mag_avg,Ax11_mag_std];
Bx11_mag = [Bx11_mag_avg,Bx11_mag_std];


end




function [qe,qs,qa, qds] = get_efficiencies(sphere_results)
    for lda_idx = 1:size(sphere_results,1)
        for instance_idx = 1:size(sphere_results,2)
            % Grab all the coeffs for a specific lda, instance pair
            qe(lda_idx, instance_idx) = sphere_results{lda_idx, instance_idx}.qe; 
            qs(lda_idx, instance_idx) = sphere_results{lda_idx, instance_idx}.qs;             
            qa(lda_idx, instance_idx) = sphere_results{lda_idx, instance_idx}.qa;      
        end
    end
qds = qe-qs-qa;


end

function [I0_par, I0_perp, I180_par, I180_perp] = get_I(modes_perp, modes_par)

    %% GET S AT THETA = 0, 180
    %[a011,a111,a021,a121,a221,a031,a131,a231,a331]
    a1n1_perp = squeeze(modes_perp.anm1([2,4,7],:,:));
    a1n2_perp = squeeze(modes_perp.anm2([2,4,7],:,:));
    b1n1_perp = squeeze(modes_perp.bnm1([2,4,7],:,:));
    b1n2_perp = squeeze(modes_perp.bnm2([2,4,7],:,:));

    a1n1_par = squeeze(modes_par.anm1([2,4,7],:,:));
    a1n2_par = squeeze(modes_par.anm2([2,4,7],:,:));
    b1n1_par = squeeze(modes_par.bnm1([2,4,7],:,:));
    b1n2_par = squeeze(modes_par.bnm2([2,4,7],:,:));
    %% GET S FOR THETA = 180 

    S1_180 = zeros(1,size(a1n1_perp,2));
    S2_180 = zeros(1,size(a1n1_perp,2));
    S3_180 = zeros(1,size(a1n1_perp,2));
    S4_180 = zeros(1,size(a1n1_perp,2));
    for n = 1:3
        S1_180 = S1_180+((-1i).^n).*((-1).^n).*(2*n+1).*(b1n1_perp(n,:,:)-a1n2_perp(n,:,:));
        S2_180 = S2_180+((-1i).^(n+1)).*((-1).^n).*(2*n+1).*(b1n2_par(n,:,:)-a1n1_par(n,:,:));
        S3_180 = S3_180+((-1i).^(n+1)).*((-1).^n).*(2*n+1).*(b1n2_perp(n,:,:)-a1n1_perp(n,:,:));
        S4_180 = S4_180+((-1i).^n).*((-1).^n).*(2*n+1).*(b1n1_par(n,:,:)-a1n2_par(n,:,:));
    end

    %% GET S FOR THETA = 0 

    S1_0 = zeros(1,size(a1n1_perp,2));
    S2_0 = zeros(1,size(a1n1_perp,2));
    S3_0 = zeros(1,size(a1n1_perp,2));
    S4_0 = zeros(1,size(a1n1_perp,2));
    for n = 1:3
        S1_0 = S1_0+((-1i).^n).*((-1)).*(2*n+1).*(b1n1_perp(n,:,:)+a1n2_perp(n,:,:));
        S2_0 = S2_0+((-1i).^(n+1)).*((-1)).*(2*n+1).*(b1n2_par(n,:,:)+a1n1_par(n,:,:));
        S3_0 = S3_0+((-1i).^(n+1)).*((-1)).*(2*n+1).*(b1n2_perp(n,:,:)+a1n1_perp(n,:,:));
        S4_0 = S4_0+((-1i).^n).*((-1)).*(2*n+1).*(b1n1_par(n,:,:)+a1n2_par(n,:,:));
    end
    %%
    S1_0 = squeeze(S1_0);
    S2_0 = squeeze(S2_0);
    S3_0 = squeeze(S3_0);
    S4_0 = squeeze(S4_0); 
    S1_180 = squeeze(S1_180);
    S2_180 = squeeze(S2_180);
    S3_180 = squeeze(S3_180);
    S4_180 = squeeze(S4_180);      
     
%     figure, 
%     plot(mean(abs(S1_0),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(abs(S2_0),2),':','LineWidth',6)   
%     plot(mean(abs(S3_0),2),'-','LineWidth',6)   
%     plot(mean(abs(S4_0),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
%     figure, 
%     plot(mean(abs(S1_180),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(abs(S2_180),2),':','LineWidth',6)   
%     plot(mean(abs(S3_180),2),'-','LineWidth',6)   
%     plot(mean(abs(S4_180),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
%     
%     figure, 
%     plot(mean(angle(S1_0),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(angle(S2_0),2),':','LineWidth',6)   
%     plot(mean(angle(S3_0),2),'-','LineWidth',6)   
%     plot(mean(angle(S4_0),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
    
    for idx = 1:size(S1_0,1)
        for idx2 = 1:size(S1_0,2)
        s23_iso(idx,idx2) = abs(S2_0(idx,idx2)).^2+abs(S3_0(idx,idx2)).^2;    
        s23_0_mix(idx,idx2) = dot(S2_0(idx,idx2),conj(S3_0(idx,idx2)))+dot(conj(S2_0(idx,idx2)),(S3_0(idx,idx2)));
        s23_180_mix(idx, idx2) = dot(S2_180(idx,idx2),conj(S3_180(idx,idx2)))+dot(conj(S2_180(idx,idx2)),(S3_180(idx,idx2)));
        s14_0_mix(idx,idx2) = dot(S4_0(idx,idx2),conj(S1_0(idx,idx2)))+dot(conj(S4_0(idx,idx2)),(S1_0(idx,idx2)));
        s14_180_mix(idx,idx2) =  dot(S4_180(idx,idx2),conj(S1_180(idx,idx2)))+dot(conj(S4_180(idx,idx2)),(S1_180(idx,idx2))); 
        
        I0_perp(idx,idx2) =  abs(S2_0(idx,idx2)).^2+abs(S3_0(idx,idx2)).^2+...
            dot(S2_0(idx,idx2),conj(S3_0(idx,idx2)))+dot(conj(S2_0(idx,idx2)),(S3_0(idx,idx2)));
        
        I0_par(idx,idx2) =  abs(S4_0(idx,idx2)).^2+abs(S1_0(idx,idx2)).^2+...
            dot(S4_0(idx,idx2),conj(S1_0(idx,idx2)))+dot(conj(S4_0(idx,idx2)),(S1_0(idx,idx2)));
        
        I180_perp(idx,idx2) =  abs(S2_180(idx,idx2)).^2+abs(S3_180(idx,idx2)).^2+...
            dot(S2_180(idx,idx2),conj(S3_180(idx,idx2)))+dot(conj(S2_180(idx,idx2)),(S3_180(idx,idx2)));
        
        I180_par(idx,idx2) =  abs(S4_180(idx,idx2)).^2+abs(S1_180(idx,idx2)).^2+...
            dot(S4_180(idx,idx2),conj(S1_180(idx,idx2)))+dot(conj(S4_180(idx,idx2)),(S1_180(idx,idx2)));    
        end
    end  
    
%     figure,
%     plot(mean(s23_0_mix,2))
%     hold on 
%     plot(mean(s23_180_mix,2))
%     plot(mean(s14_0_mix,2))
%     plot(mean(s14_180_mix,2)) 
%     hold off
%     legend('S_{23}(\theta = 0)','S_{23}(\theta = 180)','S_{14}(\theta = 0)','S_{14}(\theta = 180)')
   %keyboard;
   
    
    
end

function [I0_par, I0_perp, I180_par, I180_perp] = get_I_compact(modes)

    %% GET S AT THETA = 0, 180
    %[a011,a111,a021,a121,a221,a031,a131,a231,a331]
    Amn = squeeze(modes.Amn([2,4,7],:,:));
    Bmn = squeeze(modes.Bmn([2,4,7],:,:));
    Axmn = squeeze(modes.Axmn([2,4,7],:,:));
    Bxmn = squeeze(modes.Bxmn([2,4,7],:,:));

    %% GET S FOR THETA = 180 

    S1_180 = zeros(1,size(Amn,2));
    S2_180 = zeros(1,size(Amn,2));
    S3_180 = zeros(1,size(Amn,2));
    S4_180 = zeros(1,size(Amn,2));
    for n = 1:3
        S1_180 = S1_180+((-1i).^n).*((-1).^n).*(2*n+1).*(Bmn(n,:,:)-Amn(n,:,:));
        S2_180 = S2_180+((-1i).^(n+1)).*((-1).^n).*(2*n+1).*(Bmn(n,:,:)-Amn(n,:,:));
        S3_180 = S3_180+((-1i).^(n+1)).*((-1).^n).*(2*n+1).*(Bxmn(n,:,:)-Axmn(n,:,:));
        S4_180 = S4_180+((-1i).^n).*((-1).^n).*(2*n+1).*(Bxmn(n,:,:)-Axmn(n,:,:));
    end

    %% GET S FOR THETA = 0 

    S1_0 = zeros(1,size(Amn,2));
    S2_0 = zeros(1,size(Amn,2));
    S3_0 = zeros(1,size(Amn,2));
    S4_0 = zeros(1,size(Amn,2));
    for n = 1:3
        S1_0 = S1_0+((-1i).^n).*((-1)).*(2*n+1).*(Bmn(n,:,:)+Amn(n,:,:));
        S2_0 = S2_0+((-1i).^(n+1)).*((-1)).*(2*n+1).*(Bmn(n,:,:)+Amn(n,:,:));
        S3_0 = S3_0+((-1i).^(n+1)).*((-1)).*(2*n+1).*(Bxmn(n,:,:)+Axmn(n,:,:));
        S4_0 = S4_0+((-1i).^n).*((-1)).*(2*n+1).*(Bxmn(n,:,:)+Axmn(n,:,:));
    end
    %%
    S1_0 = squeeze(S1_0);
    S2_0 = squeeze(S2_0);
    S3_0 = squeeze(S3_0);
    S4_0 = squeeze(S4_0); 
    S1_180 = squeeze(S1_180);
    S2_180 = squeeze(S2_180);
    S3_180 = squeeze(S3_180);
    S4_180 = squeeze(S4_180);      
     
%     figure, 
%     plot(mean(abs(S1_0),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(abs(S2_0),2),':','LineWidth',6)   
%     plot(mean(abs(S3_0),2),'-','LineWidth',6)   
%     plot(mean(abs(S4_0),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
%     figure, 
%     plot(mean(abs(S1_180),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(abs(S2_180),2),':','LineWidth',6)   
%     plot(mean(abs(S3_180),2),'-','LineWidth',6)   
%     plot(mean(abs(S4_180),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
%     
%     figure, 
%     plot(mean(angle(S1_0),2),'-','LineWidth',6)   
%     hold on 
%     plot(mean(angle(S2_0),2),':','LineWidth',6)   
%     plot(mean(angle(S3_0),2),'-','LineWidth',6)   
%     plot(mean(angle(S4_0),2),':','LineWidth',6)   
%     hold off
%     legend('s1','s2','s3','s4')
%     
    
    for idx = 1:size(S1_0,1)
        for idx2 = 1:size(S1_0,2)
        s23_iso(idx,idx2) = abs(S2_0(idx,idx2)).^2+abs(S3_0(idx,idx2)).^2;    
        s23_0_mix(idx,idx2) = dot(S2_0(idx,idx2),conj(S3_0(idx,idx2)))+dot(conj(S2_0(idx,idx2)),(S3_0(idx,idx2)));
        s23_180_mix(idx, idx2) = dot(S2_180(idx,idx2),conj(S3_180(idx,idx2)))+dot(conj(S2_180(idx,idx2)),(S3_180(idx,idx2)));
        s14_0_mix(idx,idx2) = dot(S4_0(idx,idx2),conj(S1_0(idx,idx2)))+dot(conj(S4_0(idx,idx2)),(S1_0(idx,idx2)));
        s14_180_mix(idx,idx2) =  dot(S4_180(idx,idx2),conj(S1_180(idx,idx2)))+dot(conj(S4_180(idx,idx2)),(S1_180(idx,idx2))); 
        
        I0_perp(idx,idx2) =  abs(S2_0(idx,idx2)).^2+abs(S3_0(idx,idx2)).^2+...
            dot(S2_0(idx,idx2),conj(S3_0(idx,idx2)))+dot(conj(S2_0(idx,idx2)),(S3_0(idx,idx2)));
        
        I0_par(idx,idx2) =  abs(S4_0(idx,idx2)).^2+abs(S1_0(idx,idx2)).^2+...
            dot(S4_0(idx,idx2),conj(S1_0(idx,idx2)))+dot(conj(S4_0(idx,idx2)),(S1_0(idx,idx2)));
        
        I180_perp(idx,idx2) =  abs(S2_180(idx,idx2)).^2+abs(S3_180(idx,idx2)).^2+...
            dot(S2_180(idx,idx2),conj(S3_180(idx,idx2)))+dot(conj(S2_180(idx,idx2)),(S3_180(idx,idx2)));
        
        I180_par(idx,idx2) =  abs(S4_180(idx,idx2)).^2+abs(S1_180(idx,idx2)).^2+...
            dot(S4_180(idx,idx2),conj(S1_180(idx,idx2)))+dot(conj(S4_180(idx,idx2)),(S1_180(idx,idx2)));    
        end
    end  
    
%     figure,
%     plot(mean(s23_0_mix,2))
%     hold on 
%     plot(mean(s23_180_mix,2))
%     plot(mean(s14_0_mix,2))
%     plot(mean(s14_180_mix,2)) 
%     hold off
%     legend('S_{23}(\theta = 0)','S_{23}(\theta = 180)','S_{14}(\theta = 0)','S_{14}(\theta = 180)')
%    keyboard;
   
    
    
end


function [modes_perp, modes_par, modes] = grab_modes(sphere_coeffs_perp, sphere_coeffs_par)

[anm1_par, anm2_par, bnm1_par, bnm2_par,x]= convert_modes(sphere_coeffs_par);
modes_par = struct(...
    'anm1', anm1_par,...
    'bnm1', bnm1_par,...
    'anm2', anm2_par,...
    'bnm2', bnm2_par,...
    'x',x);
 
[anm1_perp, anm2_perp, bnm1_perp, bnm2_perp,x]= convert_modes(sphere_coeffs_perp);
modes_perp = struct(...
    'anm1', anm1_perp,...
    'bnm1', bnm1_perp,...
    'anm2', anm2_perp,...
    'bnm2', bnm2_perp,...
    'x',x);

% modes = struct(...
%     'Amn', cat(3, anm2_perp, anm1_par),...
%     'Bmn', cat(3, bnm1_perp, bnm2_par),...
%     'Axmn', cat(3,anm1_perp, anm2_par),...
%     'Bxmn', cat(3,bnm2_perp, bnm1_par),...
%     'x',cat(2,x,x) );

modes = struct(...
    'Bmn', cat(3, anm2_perp, bnm2_par),...
    'Amn', cat(3, bnm1_perp, anm1_par),...
    'Bxmn', cat(3,anm1_perp, bnm1_par),...
    'Axmn', cat(3,bnm2_perp, anm2_par),...
    'x',cat(2,x,x) );



end

function [anm1, anm2, bnm1, bnm2, x,...
    qext_par, qsca_par, qabs_par,...
    qext_per, qsca_per, qabs_per]= convert_modes(sphere_coeffs)

%[-22,-12,02,12,22]
%[22,12,02,-12,-22]
%n = 2
% take [n+1,end] of both arrays
% add acording to m

dummy_te = {};
dummy_tm = {};

    for lda_idx = 1:size(sphere_coeffs,1)
        for instance_idx = 1:size(sphere_coeffs,2)
            % Grab all the coeffs for a specific lda, instance pair
            dummy_te = sphere_coeffs{lda_idx, instance_idx}.a_te;
            dummy_tm = sphere_coeffs{lda_idx, instance_idx}.a_tm;   
            x(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.ka;
            qsca_par(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qsca_par; 
            qabs_par(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qabs_par;             
            qext_par(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qext_par;             
            qsca_per(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qsca_per; 
            qabs_per(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qabs_per;             
            qext_per(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qext_per; 
            
            
            
            for n = 1:length(dummy_te)
                if n == 1
                   % [-11,01,11]
                   dummy_tm_n = dummy_te{n};
                   dummy_te_n = dummy_tm{n};

                   a111 = (dummy_tm_n(3)-dummy_tm_n(1));
                   a011 = dummy_tm_n(2);

                   a112 = 2.*(dummy_te_n(3)-dummy_te_n(1));
                   a012 = dummy_te_n(2); 

                   b111 = (dummy_tm_n(3)+dummy_tm_n(1));
                   b011 = dummy_tm_n(2);

                   b112 = 2.*(dummy_te_n(3)+dummy_te_n(1));
                   b012 = dummy_te_n(2);               
                end

                if n == 2
                   %[-22,-12,02,12,22]                
                   dummy_tm_n = dummy_te{n};
                   dummy_te_n = dummy_tm{n};

                   a221 = (dummy_tm_n(5)+dummy_tm_n(1));
                   a121 = (dummy_tm_n(4)-dummy_tm_n(2));              
                   a021 = dummy_tm_n(3);

                   a222 = 2.*(dummy_te_n(5)+dummy_te_n(1));
                   a122 = 2.*(dummy_te_n(4)-dummy_te_n(2));              
                   a022 = dummy_te_n(3);               


                   b221 = (dummy_tm_n(5)-dummy_tm_n(1));
                   b121 = (dummy_tm_n(4)+dummy_tm_n(2));              
                   b021 = dummy_tm_n(3);

                   b222 = 2.*(dummy_te_n(5)-dummy_te_n(1));
                   b122 = 2.*(dummy_te_n(4)+dummy_te_n(2));              
                   b022 = dummy_te_n(3);              
                end            

                if n == 3
                   %[-33,-23,-13,03,13,23,33]                
                   dummy_tm_n = dummy_te{n};
                   dummy_te_n = dummy_tm{n};

                   a331 = (dummy_tm_n(7)-dummy_tm_n(1));
                   a231 = (dummy_tm_n(6)+dummy_tm_n(2));
                   a131 = (dummy_tm_n(5)-dummy_tm_n(1));              
                   a031 = dummy_tm_n(4);

                   a332 = 2.*(dummy_te_n(7)-dummy_te_n(1));
                   a232 = 2.*(dummy_te_n(6)+dummy_te_n(2));
                   a132 = 2.*(dummy_te_n(5)-dummy_te_n(1));              
                   a032 = dummy_te_n(4);              

                   b331 = (dummy_tm_n(7)+dummy_tm_n(1));
                   b231 = (dummy_tm_n(6)-dummy_tm_n(2));
                   b131 = (dummy_tm_n(5)+dummy_tm_n(1));              
                   b031 = dummy_tm_n(4);

                   b332 = 2.*(dummy_te_n(7)+dummy_te_n(1));
                   b232 = 2.*(dummy_te_n(6)-dummy_te_n(2));
                   b132 = 2.*(dummy_te_n(5)+dummy_te_n(1));              
                   b032 = dummy_te_n(4);                             
                end               
            end

            anm1(:,lda_idx,instance_idx) = [a011,a111,a021,a121,a221,a031,a131,a231,a331]; 
            anm2(:,lda_idx,instance_idx) = [a012,a112,a022,a122,a222,a032,a132,a232,a332]; 
            bnm1(:,lda_idx,instance_idx) = [b011,b111,b021,b121,b221,b031,b131,b231,b331]; 
            bnm2(:,lda_idx,instance_idx) = [b012,b112,b022,b122,b222,b032,b132,b232,b332];         

        end
    end

end
