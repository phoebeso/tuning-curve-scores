% Determines the head directions cells, which are defined as the cells in
% the 99th percentile of mean vector length of the spikerates per 6 degree head
% direction bin 

% bin to 0.5 s instead of 33 ms, make time bins bigger for more accuracy? 

clear all; clc; 

files = dir('SargoliniMoser2006');

nHdBins = 60;

% col 1 is cell file name, col 2 is direction, col 3 is spiketrain, col 4 is hd score
cellData = cell(length(files), 4);

for nFile = 1:length(files)
    file = files(nFile);
    filename = file.name;
    [~,name,ext] = fileparts(filename);
    if (~strcmp(ext, '.mat'))
        continue
    end    
    
    fullFileName = fullfile('SargoliniMoser2006', filename);
    load(fullFileName)
    
    dt = t(3)-t(2);
    timebins = [t; (t(end) + dt)];
    spiketrain = histcounts(ts, timebins)'; 
    
    direction = atan2(y2-y1,x2-x1)+pi/2;
    direction(direction < 0) = direction(direction<0)+2*pi; % go from 0 to 2*pi, without any negative numbers
    
    % calculate hd tuning curve and score
    [hdOccupancy, hdRates, hdScore] = calculate_hd_score(direction,spiketrain,dt,nHdBins); % hd score
    
    % Store cell data for shuffling procedure
    cellData{nFile, 1} = name; 
    cellData{nFile, 2} = direction;
    cellData{nFile, 3} = spiketrain;
    cellData{nFile, 4} = hdScore;
    
    % Graph and save tuning and occupancy curves
    hdBins = (0:2*pi/nHdBins:2*pi)';
    
    figure1 = figure(1);
    subplot('Position', [0.05 0.30 0.40 0.40])
    hdCurve = [hdRates; hdRates(1)];
    polar(hdBins, hdCurve);
    xlabel(sprintf('Head Direction Score: %f', hdScore))
    title('Head Direction Tuning Curve')
    
    subplot('Position', [0.55 0.30 0.40 0.40])
    hdCount = [hdOccupancy; hdOccupancy(1)];
    polar(hdBins, hdCount);
    title('Amount of time rat raced each direction')
    
    mkdir(['Sargolini Output/' name])
    saveas(figure1,[pwd sprintf('/Sargolini Output/%s/Head Direction 1.fig',name)]);
    close all
    
end

% Removes empty cells from cell array 
cellData(all(cellfun(@isempty,cellData),2), : ) = [];
shiftedHdScores = zeros(length(cellData)*100,1);

% Shift and calculate hd score for each cell 100 times
for i = 1:length(cellData)
    direction2 = cellData{i, 2};
    spiketrain2 = cellData{i, 3};
    
    for j = 1:100
        minShift = ceil(20/dt); % min shift is 20 s
        maxShift = length(spiketrain2)-(20/dt); % max shift is length of trial minus 20 s 
        randShift = round(minShift + rand*(maxShift-minShift)); % amount to shift spiketrain by
        
        shiftedSpiketrain = circshift(spiketrain2,randShift); % shifted spiketrain
        
        idx = ((i-1)*100)+j;
        [~, ~, shiftedScore] = calculate_hd_score(direction2,shiftedSpiketrain,dt,nHdBins); % shifted hd score
        shiftedHdScores(idx) = shiftedScore; 
    end
    
end


% determines 95th percentile of cells 
sigPercentile = prctile(shiftedHdScores,95);
hdCellsIdx = find(cell2mat(cellData(:,4)) > sigPercentile);
hdCells = cellData(hdCellsIdx, 1);

save('Sargolini Output/hd_cell_data.mat', 'hdCells', 'sigPercentile', 'cellData')
