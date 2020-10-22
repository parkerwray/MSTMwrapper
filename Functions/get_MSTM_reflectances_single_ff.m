function total_reflectances = get_MSTM_reflectances_single_ff(I, qsistats, qsdstats, ...
    center_radiis, wopt, fill_fraction)
% Formula:ff*(qsi+qsd) / (1 + FBR)
%sca = 

size_I = size(I);
if (size_I(2) ~= 1 || size_I(3) ~= 1)
   throw(MException ("Inputs:InvalidSizeError", ...
       "Invalid size of far-field intensity array"));
end
if length(fill_fraction) ~= 1
    throw(MException("Inputs:InvalidSizeError", ...
        "Invalid size of fill fraction. Must be single value"));
end
if length(wopt) ~= length(center_radiis) || length(wopt) ~= size(I, 1)
    throw(MException("Inputs:InvalidSizeError", ...
        "Number of radii is different between variables"));
end

for i = 1:length(center_radiis)
    sca(i, :) = fill_fraction * wopt(i) * ...
        (squeeze(qsistats.mean(i,1, :) + qsdstats.mean(i, 1, :)));
    fbr_num(i, :) = wopt(i) * I(i, 1, :, :, 1);
    fbr_denom(i, :) = wopt(i) * I(i, 1, :, :, 2);
end
fbr_net = sum(fbr_num, 1) ./ sum(fbr_denom, 1);
sca_net = sum(sca, 1);
total_reflectances = squeeze(sca_net ./ (1 + fbr_net));

end