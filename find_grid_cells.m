clear all; clc

% load data_for_cell77.mat
% 
% posx = posx_c;
% posy = posy_c;
% nPosBins = 20;
% boxSize = 100;
% 
% filter = gaussian(-4:4,2,0); filter = filter / sum(filter); 
% dt = post(3) - post(2); firingRate = spiketrain / dt;
% smoothFiringRate = conv(firingRate,filter,'same');
% 
% [rateMap] = compute_2d_tuning_curve(posx,posy,smoothFiringRate,nPosBins,0,boxSize);

load simdata.mat

for i = 1:4

    rateMap = simdata{i};

    [~, gridScore] = find_central_field(rateMap);

    matrix = calculate_correlation_matrix(rateMap);
    
    [rotations, correlations, circularMatrix] = calculate_spatial_periodicity(matrix);
    
    figure(i)
    subplot(2,2,1)
    imagesc(rateMap); colorbar
    title('Firing Rate Map')
    xlabel(['Grid Score: ' num2str(gridScore)])

    subplot(2,2,2)
    imagesc(matrix, [-1 1]); colorbar
    axis off
    title('Autocorrelation Matrix')

    subplot(2,2,3)
    imagesc(circularMatrix, [-1 1]); colorbar
    axis off
    title('Circular Autocorrelation Matrix')

	subplot(2,2,4)
    plot(rotations, correlations);
    xlim([0 360])
    xlabel('Rotation (deg)')
    ylabel('Correlation')
    title('Periodicity')

end
