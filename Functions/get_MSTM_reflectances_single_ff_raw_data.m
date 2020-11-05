function total_reflectances = get_MSTM_reflectances_single_ff_raw_data(...
    Ipar, Iper, qsi, qsd, center_radiis, wopt, fill_fraction)
% I is NOT squeezed, nor are the qis
% Formula:ff*(qsi+qsd) / (1 + FBR)
%sca = 

I = (mean(Ipar,3)+mean(Iper,3))./2;

% NOTE: DIMENSIONS ARE HARDCODED HERE!!!!
% This function in general is only compatible with the one format so it
% should be fine
[qsistats] = get_statistics(qsi, 3);
[qsdstats] = get_statistics(qsd, 3);


total_reflectances = get_MSTM_reflectances_single_ff(I, qsistats, qsdstats, ...
    center_radiis, wopt, fill_fraction);

end