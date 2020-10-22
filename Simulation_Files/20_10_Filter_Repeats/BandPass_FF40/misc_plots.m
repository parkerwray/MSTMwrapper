dim_Ndistributions = 2;
dim_wavelengths = 3;

[Astats] = get_mode_statistics(A, dim_Ndistributions);
[Bstats] = get_mode_statistics(B, dim_Ndistributions);
[Axstats] = get_mode_statistics(Ax, dim_Ndistributions);
[Bxstats] = get_mode_statistics(Bx, dim_Ndistributions);

% Grab the range min/max of the average mode energy in the wavelength
% dimension. We can use this to determine if a mode is negligable. 
Astats.energy = get_range(Astats.mean.mag,dim_wavelengths);
Bstats.energy = get_range(Bstats.mean.mag,dim_wavelengths);
Axstats.energy = get_range(Axstats.mean.mag,dim_wavelengths);
Bxstats.energy = get_range(Bxstats.mean.mag,dim_wavelengths);











%%
I = squeeze((mean(Ipar,3)+mean(Iper,3))./2);
%fbr = (I(:, :, :, 1)./I(:, :, :, end));
fbr = (I(:, :, :, 1)./I(:, :, :, end));
%%
I2 = squeeze((mean(Ipar(:,:,5:15,:,:),3)+mean(Iper(:,:,5:15,:,:),3))./2);
%fbr = (I(:, :, :, 1)./I(:, :, :, end));
fbr2 = (I2(:, :, 1)./I2(:, :, end));
%%
figure,
hold on
plot(wavelengths, squeeze(log10(fbr(1:4,4,:))));
%plot(wavelengths, squeeze(log10(fbr2(4,:))));
%legend("full", "partial");
%%
%plot(wavelengths, squeeze(log10(fbr(1,:,1))));
dim_N = 3;
x2 = (2*pi*44 ./ wavelengths).^2.';
scaled_A = 2 * 3 * squeeze(mean(abs(A(4,:,:,1,2)).^2,dim_N)) ./ x2;
scaled_B = 2 * 3 * squeeze(mean(abs(B(4,:,:,1,2)).^2,dim_N)) ./ x2;
%scaled_B = 2 * 3 * squeeze(mean(B(:,:,1,2),1)) ./ x2;
figure,
hold on
plot(wavelengths, squeeze(scaled_B));
plot(wavelengths, squeeze(scaled_A));
%xlim([375,500]);
hold off
%%

x2 = ((2*pi*44 ./ wavelengths).^2).';
x2 = repmat(x2, 1,5).';
scaled_A = 3 * squeeze(A(1:5,:,1,2)) ./ x2;
scaled_B = 3 * squeeze(B(1:5,:,1,2)) ./ x2;
figure,
hold on
plot(wavelengths, squeeze(abs(scaled_B)));
plot(wavelengths, squeeze(abs(scaled_A)));
%xlim([375,500]);
hold off

%%
idx_cr = 3;
idx_ff = 4;
figure,
hold on
plot(wavelengths, 6 * squeeze(Astats.mean.mag(idx_cr, idx_ff, 1,:,1,2)) ./x2 );
plot(wavelengths, 6 * squeeze(Bstats.mean.mag(idx_cr, idx_ff, 1,:,1,2)) ./x2 );
plot(wavelengths, 6 * squeeze(Axstats.mean.mag(idx_cr, idx_ff, 1,:,1,2)) ./x2 );
plot(wavelengths, 6 * squeeze(Bxstats.mean.mag(idx_cr, idx_ff, 1,:,1,2)) ./x2 );
%xlim([375,500]);
%%
figure,
hold on
shadedErrorBar(wavelengths, 6 * squeeze(Astats.mean.mag(4,1,:,1,2)) ./ x2, ...
    6 * squeeze(Astats.std.mag(4,1,:,1,2)) ./ x2)
shadedErrorBar(wavelengths, 6 * squeeze(Bstats.mean.mag(4,1,:,1,2)) ./ x2, ...
    6 * squeeze(Bstats.std.mag(4,1,:,1,2)) ./ x2)
%xlim([375, 500]);

%%

figure,
hold on
idx_r = 6;
idx_ff = 4;
shadedErrorBar(wavelengths, squeeze(Astats.mean.phase(idx_ff,1,:,1,2)), ...
    squeeze(Astats.std.phase(idx_ff,1,:,1,2)))
shadedErrorBar(wavelengths, squeeze(Bstats.mean.phase(idx_ff,1,:,1,2)), ...
    squeeze(Bstats.std.phase(idx_ff,1,:,1,2)))
shadedErrorBar(wavelengths, squeeze(Axstats.mean.phase(idx_ff,1,:,1,2)), ...
    squeeze(Axstats.std.phase(idx_ff,1,:,1,2)))
shadedErrorBar(wavelengths, squeeze(Bxstats.mean.phase(idx_ff,1,:,1,2)), ...
    squeeze(Bxstats.std.phase(idx_ff,1,:,1,2)))
%xlim([375, 500]);
%%
figure,
plot(wavelengths, squeeze(real(A(4,1:2,:,1,2))));
%% Calculate Reflectance
% Formula:ff*(qsi+qsd) / (1 + FBR)
%sca = 
I = squeeze((mean(Ipar,3)+mean(Iper,3))./2);
reflectances = zeros(length(center_radiis), length(fill_fractions), length(wavelengths));
for i = 1:length(center_radiis)
    for j = 1:length(fill_fractions)
        sca(i, j, :) = fill_fractions(j) * wopt(i) * ...
            (squeeze(qsistats.mean(i, j, 1, :) + qsdstats.mean(i, j, 1, :)));
        fbr_num(i, j, :) = wopt(i) * I(i, j, :, 1);
        fbr_denom(i, j, :) = wopt(i) * I(i, j, :, 2);
    end
end
fbr_net = sum(fbr_num, 1) ./ sum(fbr_denom, 1);
sca_net = sum(sca, 1);
total_reflectances = squeeze(sca_net ./ (1 + fbr_net));

%%
figure,
hold on
plot(wavelengths, 100*squeeze(total_reflectances(4, :)), 'LineWidth', 6);
%plot(wavelengths, squeeze(sca_net))
%plot(wavelengths, squeeze(qastats.mean(1:4,4,1,:)));
plot(lda, 100*Ropt, 'LineWidth', 6);
plot(lda, 100*Rideal, 'LineWidth', 6);
xlabel('Wavelength');
ylabel('Reflectance (%)');
ylim([0, 100]);
title('40% Fill Fraction Band-Pass Filter');
legend("True data", "Isolated data", "Ideal");
set(gca, 'FontSize', 24);
       box on
       pbaspect([1,1,1]);
       %ylim([0, 100]);