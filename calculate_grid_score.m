function [gridScore] = calculate_grid_score(rateMap)
% Calculates the grid score of a cell by rotating the rate map of the cell
% and computing the correlation between the rotated map and the original

% Expected peak correlations of sinusoidal modulation
correlation60 = calculate_correlation(rateMap, 60);
correlation120 = calculate_correlation(rateMap, 120);

% Expected trough correlations 
correlation30 = calculate_correlation(rateMap, 30);
correlation90 = calculate_correlation(rateMap, 90);
correlation150 = calculate_correlation(rateMap, 150);

% 'Gridness' of cell defined as difference between lowest expected peak
% correlation and highest expected trough correlation 
peaks = [correlation60 correlation120];
troughs = [correlation30 correlation90 correlation150];
gridScore = min(peaks) - max(troughs);

end
