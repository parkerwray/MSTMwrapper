

core = '/home/elipaul';
addpath(genpath(strcat(core,'/hypnos/Codes/MSTMwrapper')));
addpath(genpath(strcat(core,'/hypnos/Codes/Matlab_Functions')));
reflectance_func = @(A,Ax,B,Bx,Ipar,Iper,qa, qe, qsi, qsd, radii, wopt, ff) ...
    get_MSTM_reflectances_single_ff_raw_data(Ipar, Iper, qsi, qsd, radii, wopt, ff);
%%
N=1000;

[m1, s1] = bootstrap_data(N, 50, 5, A, Ax, B, Bx, Ipar, Iper, qa, qe, qsi, qsd, ...
    center_radiis, wopt, fill_fractions, reflectance_func);

[m2, s2] = bootstrap_data(N, 50, 15, A, Ax, B, Bx, Ipar, Iper, qa, qe, qsi, qsd, ...
    center_radiis, wopt, fill_fractions, reflectance_func);

[m3, s3] = bootstrap_data(N, 50, 30, A, Ax, B, Bx, Ipar, Iper, qa, qe, qsi, qsd, ...
    center_radiis, wopt, fill_fractions, reflectance_func);
%%
[m4, s4] = bootstrap_data(2, 50, 50, A, Ax, B, Bx, Ipar, Iper, qa, qe, qsi, qsd, ...
    center_radiis, wopt, fill_fractions, reflectance_func);

%%
figure,
hold on
shadedErrorBar(wavelengths/1000, m1*100, s1*100, 'LineProps', {'-r', 'LineWidth', 6})
shadedErrorBar(wavelengths/1000, m2*100, s2*100, 'LineProps', {'-g', 'LineWidth', 6})
shadedErrorBar(wavelengths/1000, m3*100, s3*100, 'LineProps', {'-b', 'LineWidth', 6})
plot(wavelengths/1000, m4*100, ":k", 'Linewidth', 6)
hold off
set(gca, 'FontSize', 24);
box on
pbaspect([1,1,1]);
title("\fontsize{18}Short-Pass Filter with 40% fill\newlinefraction and medium quality \newlineconvergence settings")
leg = legend("N=5", "N=15", "N=30", "N=50");
title(leg, "\fontsize{14}Number of \newlinedistributions\newlinein a bootstrap\newlinesample")
xlabel("Wavelength (um)");
ylabel("Reflectance (%)");
legend boxoff

