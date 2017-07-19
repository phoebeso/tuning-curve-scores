function [grid_score] = calculate_grid_score(graph)
% Expected peak correlations of sinusoidal modulation
corr_one = calculate_crosscorrelation(graph, 60);
corr_two = calculate_crosscorrelation(graph, 120);

% Expected trough correlations 
corr_three = calculate_crosscorrelation(graph, 30);
corr_four = calculate_crosscorrelation(graph, 90);
corr_five = calculate_crosscorrelation(graph, 150);

% 'Gridness' of cell defined as difference between lowest expected peak
% correlation and highest expected trough correlation 
peaks = [corr_one corr_two];
troughs = [corr_three corr_four corr_five];
grid_score = min(peaks) - max(troughs);

end
