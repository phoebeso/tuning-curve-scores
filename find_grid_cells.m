clear all; clc

load data_for_cell77.mat

posx = posx_c;
posy = posy_c;
nPosBins = 20;
boxSize = 100;

filter = gaussian(-4:4,2,0); filter = filter / sum(filter); 
dt = post(3) - post(2); firingRate = spiketrain / dt;
smoothFiringRate = conv(firingRate,filter,'same');

[rateMap] = compute_2d_tuning_curve(posx,posy,smoothFiringRate,nPosBins,0,boxSize);

load simdata.mat

for i = 1:1
    rateMap = simdata{4};

    % Calculates grid score from the rate map
    [~, gridScore] = find_central_field(rateMap);

    % Calculates the correlation matrix from the rate map
    correlationMatrix = calculate_correlation_matrix(rateMap);
    
    % Determines the spatial periodicity of the correlation matrix by
    % calculating the correlation of the matrix in intervals of 6 degrees
    [rotations, correlations, circularMatrix] = calculate_spatial_periodicity(correlationMatrix);
    
    % Calculates grid score from the correlation matrix
    gridScore2 = calculate_grid_score(circularMatrix);
    
    figure(i)
    subplot(2,2,1)
    imagesc(rateMap); colorbar
    title('Firing Rate Map')
    xlabel(['Grid Score: ' num2str(gridScore)])

    subplot(2,2,2)
    imagesc(correlationMatrix, [-1 1]); colorbar
    axis off
    title('Autocorrelation Matrix')

    subplot(2,2,3)
    imagesc(circularMatrix, [-1 1]); colorbar
    axis off
    title('Circular Autocorrelation Matrix')

	subplot(2,2,4)
    plot(rotations, correlations);
    xlim([0 360])
    ylim([-inf 1])
    xlabel({'Rotation (deg)'; sprintf('Grid Score: %f', gridScore2)})
    ylabel('Correlation')
    title('Periodicity')
    
    % Partitions the correlation periodicity curve and collapses the
    % partitioned data 
    collapsePartitionData = analyze_periodicity(rotations, correlations);
    
    maxCollapseValues = zeros(length(collapsePartitionData), 1);
    figure(i*2)
    for j = 1:length(collapsePartitionData)
        subplot(3,4,j)
        nPartitions = collapsePartitionData{j,1};
        partitionData = collapsePartitionData{j,3};
        [~,~] = graph_polar_data(nPartitions, partitionData);
        
        maxCollapseValues(j) = collapsePartitionData{j,2};
    end
    
    subplot(3,4,[9 10 11 12])
    bar([3 4 5 6 7 8 9 10], maxCollapseValues);
    xlabel('Number of Partitions')
    ylabel('Max Collapsed Data Value')
    
    % Performs FFT on the correlation matrix 
    [shiftedSpectrogram, polarSpectrogram, nComponents] = fourier_transform(correlationMatrix);
    
%     figure(i*3)
%     subplot(1,2,1)
%     imagesc(shiftedSpectrogram)
%     axis off
%     title('Fourier Spectogram')
%     
%     subplot(1,2,2)
%     imagesc(polarSpectrogram)
%     axis off 
%     title('Polar Fourier Spectogram')
    
end
