% Calculates and analyzes the gridness of each cell 

clear all; clc

% Loop through and load files from folder
files = dir('SargoliniMoser2006');

% col 1 is cell file name, col 2 is smoothed rate map, col 3 is correlation
% matrix, col 4 is circular matrix, cell 5 rotational correlations, cell 6
% is grid score, cell 7 is number of periods of max collapsed partition
% value, col 8 is fourier spectrogram, col 9 is the polar spectrogram, col
% 10 is array of max power of 100 shifted cells, col 11 is max power 
cellData = cell(length(files), 11);

for nFile = 1:length(files)
    file = files(nFile);
    filename = file.name;
    [~,name,ext] = fileparts(filename);
    if (~strcmp(ext, '.mat'))
        continue
    end 
    
    % Loads data from file
    fullFileName = fullfile('SargoliniMoser2006', filename);
    load (fullFileName)

    % Adjusts position vectors so that pos values range from 0-100 instead of -50-50
    x1 = x1 + 50; x2 = x2 + 50;
    y1 = y1 + 50; y2 = y2 + 50;
    posx = (x1 + x2) ./ 2;
    posy = (y1 + y2) ./ 2;

    % Calculates the spiketrain
    dt = t(3) - t(2);
    timebins = [t; (t(end) + dt)];
    spiketrain = histcounts(ts, timebins)'; 

    boxSize = 100; % Length in cm of one side of the environment
    nPosBins = 20; % Each position bin is 5 x 5 cm 

    % Smooths firing rate
    filter = gaussian(-4:4,2,0); filter = filter / sum(filter); 
    firingRate = spiketrain / dt;
    smoothFiringRate = conv(firingRate,filter,'same');

    % Calculates the unsmoothed and smoothed rate map 
    [rateMap, smoothRateMap] = calculate_2d_tuning_curve(posx,posy,smoothFiringRate,nPosBins,0,boxSize);
    
    % Converts all NaNs to zeros otherwise future calculations get thrown
    % off, not sure if this should be done 
    rateMap(isnan(rateMap)) = 0;
    smoothRateMap(isnan(smoothRateMap)) = 0;
    
    % Calculates the correlation matrix from the rate map
    correlationMatrix = calculate_correlation_matrix(smoothRateMap);
    
    % Determines the spatial periodicity of the correlation matrix by
    % calculating the correlation of the matrix in intervals of 6 degrees
    [rotations, correlations, circularMatrix, threshold] = calculate_spatial_periodicity(correlationMatrix);
    
    % Calculates grid score from the correlation matrix
    gridScore = calculate_grid_score(circularMatrix);
    
    % Plots rate map, autocorrelation matrix, spatial periodicity, and grid
    % score data 
    figure1 = figure(1);
    subplot(2,2,1)
    imagesc(smoothRateMap); colorbar
    axis off
    title('Firing Rate Map')

    subplot(2,2,2)
    imagesc(correlationMatrix, [-1 1]); colorbar
    axis off
    title('Autocorrelation Matrix')

    subplot(2,2,3)
    imagesc(circularMatrix, [-1 1]); colorbar
    set(gca,'xtick',[],'ytick',[])
    xlabel(sprintf('Threshold: %f', threshold))
    title('Circular Autocorrelation Matrix')

	subplot(2,2,4)
    plot(rotations, correlations);
    xlim([0 360])
    ylim([-inf 1])
    xlabel({'Rotation (deg)'; sprintf('Autocorrelation Matrix Grid Score: %f', gridScore)})
    ylabel('Correlation')
    title('Periodicity')
    
    % Partitions the correlation periodicity curve into 3-10 periods and 
    % collapses/sums the partitioned data 
    collapsePartitionData = analyze_periodicity(rotations, correlations);
    
    % Plots polar histograms of the collapsed partitioned data 
    maxCollapseValues = zeros(length(collapsePartitionData), 1);
    figure2 = figure(2);
    for k = 1:length(collapsePartitionData)
        subplot(3,4,k)
        nPartitions = collapsePartitionData{k,1};
        partitionData = collapsePartitionData{k,3};
%         plot_adjusted_polar_graph % Plots graph
        [~,~] = adjusted_polar_graph(nPartitions, partitionData); % Plots partition graphs
        
        maxCollapseValues(k) = collapsePartitionData{k,2};
    end
    
    % Plots a bar graph of the maximum value of the collapsed partitioned
    % data. The correct number of periods in the spatial periodicity curve
    % should yield the highest maximum value. 
    subplot(3,4,[9 10 11 12])
    hold on
    for nBar = 1:length(maxCollapseValues)
        h = bar(nBar + 2,maxCollapseValues(nBar));
        if maxCollapseValues(nBar) ~= max(maxCollapseValues)
            set(h,'FaceColor','b');
        else
            set(h,'FaceColor','r');
            maxNumPartitions = nBar + 2; 
        end
    end
    hold off
    xlabel('Number of Partitions')
    ylabel('Max Collapsed Data Value')
    
    % Calculates the two-dimensional Fourier spectrogram 
    adjustedRateMap = rateMap - nanmean(rateMap(:));
    meanFr = sum(spiketrain) / t(end);
    maxRate = max(rateMap(:));
    maxAdjustedRate = max(adjustedRateMap(:));
    [fourierSpectrogram, polarSpectrogram, rhoMeanPower, superimpose, beforeMaxPower, afterMaxPower, nComponents, shiftedPowers] = fourier_transform(adjustedRateMap, meanFr, spiketrain, dt, posx, posy);
    
    figure3 = figure(3);
    subplot(2,2,1)
    imagesc(fourierSpectrogram)
    set(gca, 'XTick', [], 'YTick', []);
    xlabel({sprintf('Max Rate of Unsmoothed Rate Map: %.3f | Mean Rate: %.3f | Max Rate of Unsmooted Adjusted Rate Map: %.3f', maxRate, nanmean(rateMap(:)), maxAdjustedRate); ...
            sprintf('Before Normalize Max Power: %.3f | Firing Rate Mean: %.3f | After Normalize Max Power: %.3f', beforeMaxPower, meanFr, afterMaxPower)})
    title('Fourier Spectogram')
    
    subplot(2,2,2)
    imagesc(polarSpectrogram)
    set(gca, 'XTick', [], 'YTick', []);
    xlabel(['Number of Main Components: ' num2str(nComponents)])
    title('Polar Fourier Spectogram')
    
    subplot(2,2,3)
    polar(deg2rad((0:360)'),rhoMeanPower)
    title('Fourier Polar Plot')
    hold on
    polar(deg2rad(superimpose{1}), superimpose{2})
    hold off
    
    cellData{nFile, 1} = name;
    cellData{nFile, 2} = smoothRateMap;
    cellData{nFile, 3} = correlationMatrix;
    cellData{nFile, 4} = circularMatrix; 
    cellData{nFile, 5} = correlations; 
    cellData{nFile, 6} = gridScore;
    cellData{nFile, 7} = maxNumPartitions;
    cellData{nFile, 8} = fourierSpectrogram;
    cellData{nFile, 9} = polarSpectrogram;
    cellData{nFile, 10} = shiftedPowers;
    cellData{nFile, 11} = afterMaxPower; 

%     Saves and closes figures
    mkdir(['Sargolini Output/' name])
    saveas(figure1,[pwd sprintf('/Sargolini Output/%s/Grid-Grid Scoring.fig',name)]);
    saveas(figure2,[pwd sprintf('/Sargolini Output/%s/Grid-Partitions.fig',name)]);
    saveas(figure3,[pwd sprintf('/Sargolini Output/%s/Grid-Fourier Transform.fig',name)]);
    close all
end

% Removes empty cells from cell array
cellData(all(cellfun(@isempty,cellData),2), : ) = [];

gridCellsIdx = find(cell2mat(cellData(:,6)) > 0);
gridCells = cellData(gridCellsIdx, 1);

sigPercentile = prctile(cell2mat(cellData(:,10)), 95);
periodicCellsIdx = find(cell2mat(cellData(:,11)) > sigPercentile);
periodicCells = cellData(periodicCellsIdx, 1);

save('Sargolini Output/grid_cell_data.mat', 'cellData', 'gridCells', 'periodicCells')
