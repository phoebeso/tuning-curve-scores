% Calculates and analyzes the gridness of each cell 

clear all; clc

% Loop through and load files from folder
files = dir('SargoliniMoser2006');
for nFile = 1:length(files)
    file = files(nFile);
    filename = file.name;
    [~,name,ext] = fileparts(filename);
    if (~strcmp(ext, '.mat'))
        continue
    end
    
    fullFileName = fullfile('SargoliniMoser2006', filename);
    load (fullFileName)
    
    % Calculates the spiketrain, firing rate, and rate map
    x1 = x1 + 50; x2 = x2 + 50; % Adjusts position vectors
    y1 = y1 + 50; y2 = y2 + 50;
    posx_c = (x1 + x2) ./ 2;
    posy_c = (y1 + y2) ./ 2;

    timebins = [t; (t(end) + t(3) - t(2))];
    spiketrain = histcounts(ts, timebins)'; 

    boxSize = 100; % Length in cm of one side of the environment
    nPosBins = 20; % Each position bin is 5 x 5 cm 

    filter = gaussian(-4:4,2,0); filter = filter / sum(filter); 
    dt = t(3) - t(2); firingRate = spiketrain / dt;
    smoothFiringRate = conv(firingRate,filter,'same');

    [rateMap] = compute_2d_tuning_curve(posx_c,posy_c,smoothFiringRate,nPosBins,0,boxSize);
    
    rateMap(isnan(rateMap)) = 0;

    % Calculates grid score from the rate map
    [~, gridScore] = find_central_field(rateMap);

    % Calculates the correlation matrix from the rate map
    correlationMatrix = calculate_correlation_matrix(rateMap);
    
    % Determines the spatial periodicity of the correlation matrix by
    % calculating the correlation of the matrix in intervals of 6 degrees
    [rotations, correlations, circularMatrix, threshold] = calculate_spatial_periodicity(correlationMatrix);
    
    % Calculates grid score from the correlation matrix
    gridScore2 = calculate_grid_score(circularMatrix);
    
    figure1 = figure(1);
    subplot(2,2,1)
    imagesc(rateMap); colorbar
    title('Firing Rate Map')
    xlabel(['Rate Map Grid Score: ' num2str(gridScore)])

    subplot(2,2,2)
    imagesc(correlationMatrix, [-1 1]); colorbar
    axis off
    title('Autocorrelation Matrix')

    subplot(2,2,3)
    imagesc(circularMatrix, [-1 1]); colorbar
    axis off
    xlabel(sprintf('Threshold: %f', threshold))
    title('Circular Autocorrelation Matrix')

	subplot(2,2,4)
    plot(rotations, correlations);
    xlim([0 360])
    ylim([-inf 1])
    xlabel({'Rotation (deg)'; sprintf('Autocorrelation Matrix Grid Score: %f', gridScore2)})
    ylabel('Correlation')
    title('Periodicity')
    
    % Partitions the correlation periodicity curve and collapses the
    % partitioned data 
    collapsePartitionData = analyze_periodicity(rotations, correlations);
    
    maxCollapseValues = zeros(length(collapsePartitionData), 1);
    figure2 = figure(2);
    for j = 1:length(collapsePartitionData)
        subplot(3,4,j)
        nPartitions = collapsePartitionData{j,1};
        partitionData = collapsePartitionData{j,3};
        [~,~] = adjusted_polar_graph(nPartitions, partitionData); % Plots graph
        
        maxCollapseValues(j) = collapsePartitionData{j,2};
    end
    
    subplot(3,4,[9 10 11 12])
    bar([3 4 5 6 7 8 9 10], maxCollapseValues);
    xlabel('Number of Partitions')
    ylabel('Max Collapsed Data Value')
    
    % Performs FFT on the correlation matrix 
    [shiftedSpectrogram, polarSpectrogram, nComponents] = fourier_transform(correlationMatrix);
    
    figure3 = figure(3);
    subplot('Position', [0.05 0.30 0.40 0.40])
    imagesc(shiftedSpectrogram)
    axis off
    title('Fourier Spectogram')
    
    subplot('Position', [0.55 0.30 0.40 0.40])
    imagesc(polarSpectrogram)
    axis off 
    title('Polar Fourier Spectogram')
    
    % Saves and closes figures
    mkdir(['Figures/' name])
    saveas(figure1,[pwd sprintf('/Figures/%s/Grid 1.fig',name)]);
    saveas(figure2,[pwd sprintf('/Figures/%s/Grid 2.fig',name)]);
    saveas(figure3,[pwd sprintf('/Figures/%s/Grid 3.fig',name)]);
    close all
end
