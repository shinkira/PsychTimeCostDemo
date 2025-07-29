function y = nanunique(x)

inds = not(isnan(x));
y = unique(x(inds));
