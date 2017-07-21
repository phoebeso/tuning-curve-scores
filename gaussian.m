function [y] = gaussian(x,sig,c)
% Calculates gaussian of inputed x values given sigma and c values
% Adapted from gaussmf function from Fuzzy Logic Toolbox 

y = exp((-(x - c).^2) / (2 * (sig^2)));

end
