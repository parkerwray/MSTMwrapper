function [y] = get_range(x, dim)

y.max = max(x,[],dim);
y.min = min(x,[],dim);

end









