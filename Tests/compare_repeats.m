%% Calculate Reflectance
% Formula:ff*(qsi+qsd) / (1 + FBR)
%sca = 
I = squeeze((mean(Ipar,3)+mean(Iper,3))./2);
reflectances = zeros(length(center_radiis), length(fill_fractions), length(wavelengths));
for i = 1:length(center_radiis)
    for j = 1:length(fill_fractions)
        sca(i, j, :) = fill_fractions(j) * wopt(i) * ...
            (squeeze(qsistats.mean(i, j, 1, :) + qsdstats.mean(i, j, 1, :)));
        fbr_num(i, j, :) = wopt(i) * I(i, :, 1);
        fbr_denom(i, j, :) = wopt(i) * I(i, :, 2);
    end
end
fbr_net = sum(fbr_num, 1) ./ sum(fbr_denom, 1);
sca_net = sum(sca, 1);
total_reflectances = squeeze(sca_net ./ (1 + fbr_net));

%%
figure,
hold on
plot(wavelengths_high, 100*squeeze(total_reflectances_high(:)), 'LineWidth', 6);
plot(wavelengths, 100*squeeze(total_reflectances(:)), 'LineWidth', 6);
%plot(wavelengths, squeeze(sca_net))
%plot(wavelengths, squeeze(qastats.mean(1:4,4,1,:)));
%plot(lda, 100*Ropt, 'LineWidth', 6);
%plot(lda, 100*Rideal, 'LineWidth', 6);
xlabel('Wavelength');
ylabel('Reflectance (%)');
ylim([0, 100]);
title('40% Fill Fraction Long-Pass Filter');
legend("High Convergence Data", "Low Convergence data");
set(gca, 'FontSize', 24);
       box on
       pbaspect([1,1,1]);
       %ylim([0, 100]);
       
       
 %%
 total_reflectances_high = total_reflectances;
 A_high  = A;
 Ax_high = Ax;
 B_high  = B;
 Bx_high = Bx;
 Astats_high  = Astats;
 Axstats_high = Axstats;
 Bstats_high  = Bstats;
 Bxstats_high = Bxstats;
 Ipar_high = Ipar;
 Iper_high = Iper;
 qa_high  = qa;
 qe_high  = qe;
 qsd_high = qsd;
 qsi_high = qsi;
 qastats_high  = qastats;
 qestats_high  = qestats;
 qsistats_high = qsistats;
 qsdstats_high = qsdstats;
 reflectances_high = reflectances;
 wavelengths_high = wavelengths;
 
 %%
 
figure,
hold on
plot(wavelengths, squeeze(abs(A(1,40,:,1,1))))
plot(wavelengths_high, squeeze(abs(A_high(1,40,:,1,1))))
legend("Original", "High Convergence")
legend boxoff

%%
figure,
hold on
plot(wavelengths, squeeze(mean(abs(A(1,:,:,1,1)),2)));
plot(wavelengths, squeeze(mean(abs(A_high(1,:,:,1,1)),2)));