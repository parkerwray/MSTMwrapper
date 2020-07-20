function [out] = get_statistics(c, dim)

%{
    This function calculates the mean and standard deviation in particle
    efficiency 
%}

out.mean =  mean(c, dim);
out.std = std(c_mag, [], dim);


end

