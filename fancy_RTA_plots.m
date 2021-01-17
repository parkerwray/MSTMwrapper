%% Calculate Reflectance
% Formula:ff*(qsi+qsd) / (1 + FBR)
%sca = 
I = squeeze((mean(Ipar,3)+mean(Iper,3))./2);
reflectances = zeros(length(center_radiis), length(fill_fractions), length(wavelengths));
for i = 1:length(center_radiis)
    sca(i, :) = fill_fractions * wopt(i) * ...
        (squeeze(qsistats.mean(i, 1, 1, :) + qsdstats.mean(i, 1, 1, :)));
    fbr_num(i, :) = wopt(i) * I(i, :, 1);
    fbr_denom(i, :) = wopt(i) * I(i, :, 2);
end
fbr_net = sum(fbr_num, 1) ./ sum(fbr_denom, 1);
sca_net = sum(sca, 1);
total_reflectances = squeeze(sca_net ./ (1 + fbr_net));
%% Calculate Absorption
% Formula: ff * qa
for i = 1:length(center_radiis)
    abs(i,:) = fill_fractions * wopt(i) * squeeze(qastats.mean(i,1,1,:));
end
total_absorptions = sum(abs,1);

%% Calculate Transmittance
% Formula: 1-reflectance-absorption
total_transmittances = 1 - total_reflectances - total_absorptions;

%% Plot fancily
um_wavelengths = wavelengths / 1000;
figure,
hold on
plot(um_wavelengths, 100*squeeze(total_reflectances(:)), '-b', 'LineWidth', 6);
plot(um_wavelengths, 100*squeeze(total_absorptions(:)), '-r', 'LineWidth', 6);
plot(um_wavelengths, 100*squeeze(total_transmittances(:)), '-k', 'LineWidth', 6);
xlabel(['Wavelength (', char(956), 'm)']);
ylabel('Reflectance, Absorbance, Transmittance (%)');
ylim([0, 100]);
thetitle = title('40% fill fraction optimized Ge long-pass filter with high convergence');
set(thetitle, "FontWeight", "normal");
legend("Reflectance", "Absorbance", "Transmittance");
legend boxoff
set(gca, 'FontSize', 24);
       box on
       pbaspect([1,1,1]);
       %ylim([0, 100]);

