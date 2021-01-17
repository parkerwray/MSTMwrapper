function [ropt_new, wopt_new] = aggregate_similar_particles(ropt, wopt, min_wopt, min_diff)

% Inputs: 
%   ropt: radii of particles
%   wopt: (area) weights of particles, i.e. the proportion of the film
%       which is occupied at all that is occupied by particles of the
%       correspongin wopt
%   min_wopt: The minimum probability to be consider important enough to
%       compute
%   min_diff: The minimum difference between two particles for them to be
%       considered different. The subtlety of the way this works is
%       described below.

% Outputs:
%   ropt_new: The new, reduced radii to consider
%   wopt_new: The corresponding weights to these radii

% To avoid overly grouping systems with high numbers of particles, we
% group particles into categories of "same" particles with span at most
% min_diff--i.e. placing 1,1.5,2.25 with min_diff of one will group [1,1.5]
% together but 2.25 separately, because 2.25-1 > 1. We sort by ascending
% values of ropt for consistency.
% We could consider recursing this function multiple times--this might help
% for the edge cases such as [1,1.9,2,2.1,2.2] --> 1.775 (instead of 
% [1.45, 2.1] as it is now) by grouping them all together when they are 
% reasonably close, but ignoring cases such as [1,2,3,4]-->[1.5,3.5] either
% way (where if just individual differences were considered it would map
% just to 2.5). This is likely a negligible difference, however.

if length(ropt) ~= length(wopt)
    throw(MException("Inputs:Mismatch", "Radii and Weight sizes do not match"));
end

[ropt, sort_idx] = sort(ropt, 'ascend');
wopt = wopt(sort_idx);

ropt_cell = {};
ropt_cell{1} = ropt(1);
wopt_cell = {};
wopt_cell{1} = wopt(1);

% Break into categories with approximately identical radii
for idx_1 = 2:length(ropt)
    never_same = true;
    for idx_2 = 1:length(ropt_cell)
        same = true;
        for idx_3 = 1:length(ropt_cell{idx_2})
            if (abs(ropt(idx_1) - ropt_cell{idx_2}(idx_3)) > min_diff)
                same = false;
            end
        end
        if same && never_same
            never_same = false;
            ropt_cell{idx_2}(idx_3 + 1) = ropt(idx_1);
            wopt_cell{idx_2}(idx_3 + 1) = wopt(idx_1);
        end
    end
    if never_same
        [ropt_cell{idx_2+1}] = ropt(idx_1);
        [wopt_cell{idx_2+1}] = wopt(idx_1);
    end
end

% Sum weights over categories and note which are negligible
non_negligible = 1:length(ropt_cell);
ropt_new = zeros(length(ropt_cell), 1);
wopt_new = zeros(length(ropt_cell), 1);
for idx = 1:length(ropt_cell)
    ropt_new(idx) = mean(ropt_cell{idx});
    wopt_new(idx) = sum (wopt_cell{idx});
    if wopt_new(idx) < min_wopt
        non_negligible(non_negligible==idx) = [];
    end
end

% Clear radii with negligible weights
ropt_new = ropt_new(non_negligible);
wopt_new = wopt_new(non_negligible);

    
end