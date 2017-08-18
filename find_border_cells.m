clear all; clc

files = dir('SargoliniMoser2006');

boxSize = 100; % Length in cm of one side of the environment
nPosBins = 20; % Each position bin is 5 x 5 cm 

% col 1 is cell file name, col 2 is posx, col 3 is posy, col 4 is
% spiketrain, col 5 is border score
cellData = cell(length(files), 5);

% Loops through files from folder
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

    % Smooths firing rate
    filter = gaussian(-4:4,2,0); filter = filter / sum(filter); 
    firingRate = spiketrain / dt;
    smoothFiringRate = conv(firingRate,filter,'same');

    % Calculates the unsmoothed and smoothed rate map
    % Smoothed rate map used for correlation calculations, unsmoothed rate
    % map used for fourier analysis
    [~, smoothRateMap] = calculate_2d_tuning_curve(posx,posy,smoothFiringRate,nPosBins,0,boxSize);
    
    smoothRateMap(isnan(smoothRateMap)) = 0; 

    % Calcultes border score
    borderScore = calculate_border_score(smoothRateMap);

    % Plots and saves data
    figure1 = figure('Name', sprintf('%s Border Score', name), 'Numbertitle', 'off');
    imagesc(smoothRateMap); colorbar
    xlabel(['Border Score: ' num2str(borderScore)])
    
    % Store data
    cellData{nFile, 1} = name;
    cellData{nFile, 2} = posx;
    cellData{nFile, 3} = posy; 
    cellData{nFile, 4} = spiketrain;
    cellData{nFile, 5} = borderScore;
    
    % Saves and closes figures
    mkdir(['Sargolini Output/' name])
    saveas(figure1,[pwd sprintf('/Sargolini Output/%s/Border Scoring.fig',name)]);
    close all
end

% Removes empty cells from cell array 
cellData(all(cellfun(@isempty,cellData),2), : ) = [];
shiftedBorderScores = zeros(length(cellData)*100,1);

% Shift and calculate border score for each cell 100 times
for i = 1:length(cellData)
    posx2 = cellData{i, 2};
    posy2 = cellData{i, 3};
    spiketrain2 = cellData{i, 4};
    
    for j = 1:100
        minShift = ceil(20/dt); % min shift is 20 s
        maxShift = length(spiketrain2)-(20/dt); % max shift is length of trial minus 20 s 
        randShift = round(minShift + rand*(maxShift-minShift)); % amount to shift spiketrain by
        
        shiftedSpiketrain = circshift(spiketrain2,randShift); % shifted spiketrain
        
        % Smooths firing rate
        filter = gaussian(-4:4,2,0); filter = filter / sum(filter); 
        shiftedFiringRate = shiftedSpiketrain / dt;
        shiftedSmoothFiringRate = conv(shiftedFiringRate,filter,'same');
        
        idx = ((i-1)*100)+j;
        [~, shiftedRateMap] = calculate_2d_tuning_curve(posx2,posy2,shiftedSmoothFiringRate,nPosBins,0,boxSize);
        shiftedRateMap(isnan(shiftedRateMap)) = 0; 
        shiftedScore = calculate_border_score(shiftedRateMap);
        shiftedBorderScores(idx) = shiftedScore; 
    end
    
end

% Determines 95th percentile of shuffled border scores to define/identify 
% border cells
sigPercentile = prctile(shiftedBorderScores,95);
borderCellsIdx = find(cell2mat(cellData(:,5)) > sigPercentile);
borderCells = cellData(borderCellsIdx, 1);

save('Sargolini Output/border_cell_data.mat', 'borderCells', 'sigPercentile', 'cellData')
