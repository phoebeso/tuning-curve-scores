% Determines the head directions cells, which are defined as the cells in
% the 99th percentile of mean vector length of the spikerates per 6 degree head
% direction bin 

% bin to 0.5 s instead of 33 ms??? make time bins bigger for more
% accuracy...? 

clear all; clc; 

load spikePos_openfield_D4.mat
load thetaEEG_JZ1_day4_nt11.mat
load JZ1posday4.mat
posStruct6 = pos{1,4}{1,6};
posStruct8 = pos{1,4}{1,8};
spikePosStruct = JZ1spikesPosD4;
thetaStruct6 = thetaEEG_JZ1_day4_epoch6_nt11;
thetaStruct8 = thetaEEG_JZ1_day4_epoch8_nt11;
numCells = length(spikePosStruct); % number of cells being analyzed

hd_scores = zeros(numCells,2);
shifted_hd_scores = zeros(numCells*100,2);
n_hd_bins = 60;

for p = 1:numCells
    epoch = spikePosStruct(p).index(2);
    if epoch == 6
        post = posStruct6.data(:,1);
        direction = posStruct6.data(:,4);
    elseif epoch == 8
        post = posStruct8.data(:,1);
        direction = posStruct8.data(:,4);
    end
    
    sampleRate = 30;
    dt = post(3)-post(2);
    
    spiketimes = spikePosStruct(p).spikes(:,1);
    timebins = linspace(post(1),post(end)+post(3)-post(2),length(post)+1);
    spiketrain = histcounts(spiketimes,timebins)';
    
    direction(direction < 0) = direction(direction<0)+2*pi; % go from 0 to 2*pi, without any negative numbers
    
    % calculate hd score
    hd_scores(p) = calculate_hd_score(direction,spiketrain,dt,n_hd_bins); % hd score
    
    % shift and calculate hd score for each cell 100 times
    for i = 1:100
        min_shift = ceil(20/dt); % min shift is 20 s
        max_shift = length(spiketrain)-(20/dt); % max shift is length of trial minus 20 s 
        rand_shift = round(min_shift + rand*(max_shift-min_shift)); % amount to shift spiketrain by
        
        spiketrain_shift = circshift(spiketrain,rand_shift); % shifted spiketrain
        
        idx = ((p-1)*100)+i;
        shifted_hd_scores(idx,1) = p; % cell num
        shifted_hd_scores(idx,2) = calculate_hd_score(direction,spiketrain_shift,dt,n_hd_bins); % shifted hd score
    end
    
end

% determines 95th percentile of cells 
sig_percentile = prctile(shifted_hd_scores(:,2),95);
% hd_cells_idx = find(shifted_hd_scores(:,2) > sig_percentile);
% hd_cells = unique(shifted_hd_scores(hd_cells_idx,1));
hd_cells = find(hd_scores > sig_percentile);

save hd_cell_data.mat

