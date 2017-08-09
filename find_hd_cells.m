% Determines the head directions cells, which are defined as the cells in
% the 99th percentile of mean vector length of the spikerates per 6 degree head
% direction bin 

% bin to 0.5 s instead of 33 ms??? make time bins bigger for more
% accuracy? 

clear all; clc; 

files = dir('SargoliniMoser2006');

numCells = 31;

hdScores = zeros(numCells,2);
shifted_hd_scores = zeros(numCells*100,2);
nHdBins = 60;

for nFile = 1:length(files)
    file = files(nFile);
    filename = file.name;
    [~,name,ext] = fileparts(filename);
    if (~strcmp(ext, '.mat'))
        continue
    end
    
    fullFileName = fullfile('SargoliniMoser2006', filename);
    load (fullFileName)
    
    dt = t(3)-t(2);
    timebins = [t; (t(end) + dt)];
    spiketrain = histcounts(ts, timebins)'; 
    
    direction = atan2(y2-y1,x2-x1)+pi/2;
    direction(direction < 0) = direction(direction<0)+2*pi; % go from 0 to 2*pi, without any negative numbers
    
    % calculate hd tuning curve and score
    [hdOccupancy, hdRates, hdScore, x, y] = calculate_hd_score(direction,spiketrain,dt,nHdBins); % hd score
    
    hdBins = (0:2*pi/nHdBins:2*pi)';
    
    figure1 = figure(1);
    subplot('Position', [0.05 0.30 0.40 0.40])
    hdCurve = [hdRates; hdRates(1)];
    polar(hdBins, hdCurve);
    xlabel(sprintf('Head Direction Score: %f', hdScore))
    title('Head Direction Tuning Curve')
    hold on
    polar(atan(y/x), sqrt(x^2 + y^2), '*')
    
    subplot('Position', [0.55 0.30 0.40 0.40])
    hdCount = [hdOccupancy; hdOccupancy(1)];
    polar(hdBins, hdCount);
    title('Amount of time rat raced each direction')
    
    mkdir(['Figures/' name])
    saveas(figure1,[pwd sprintf('/Figures/%s/Head Direction 1.fig',name)]);
    close all
    
%     % shift and calculate hd score for each cell 100 times
%     for i = 1:100
%         min_shift = ceil(20/dt); % min shift is 20 s
%         max_shift = length(spiketrain)-(20/dt); % max shift is length of trial minus 20 s 
%         rand_shift = round(min_shift + rand*(max_shift-min_shift)); % amount to shift spiketrain by
%         
%         spiketrain_shift = circshift(spiketrain,rand_shift); % shifted spiketrain
%         
%         idx = ((p-1)*100)+i;
%         shifted_hd_scores(idx,1) = p; % cell num
%         shifted_hd_scores(idx,2) = calculate_hd_score(direction,spiketrain_shift,dt,nHdBins); % shifted hd score
%     end
    
end

% % determines 95th percentile of cells 
% sig_percentile = prctile(shifted_hd_scores(:,2),95);
% hd_cells = find(hd_scores > sig_percentile);
% 
% save hd_cell_data.mat

