% Determines the speed cells, which are defined as the cells in the 99th
% percentile of Pearson correlation between spike rate and speed  

% 0.5 s bins

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

speed_scores = zeros(numCells,1);
shifted_speed_scores = zeros(numCells*100,2);

for p = 1:numCells
    epoch = spikePosStruct(p).index(2);
    if epoch == 6
        post = posStruct6.data(:,1);
        posx_c = posStruct6.data(:,2);
        posy_c = posStruct6.data(:,3);
    elseif epoch == 8
        post = posStruct8.data(:,1);
        posx_c = posStruct8.data(:,2);
        posy_c = posStruct8.data(:,3);
    end
    
    sampleRate = 30;
    dt = post(3)-post(2);
    
    spiketimes = spikePosStruct(p).spikes(:,1);
    timebins = linspace(post(1),post(end)+dt,length(post)+1);
    spiketrain = histcounts(spiketimes,timebins)';
    
    % bins = linspace(post(1),post(end)+dt,(post(end)+dt-post(1)+0.5)/0.5); % 0.5 s timebins
    
    % spiketrain for every 0.5 s instead of 1/30 s 
    rebinned_spiketrain = zeros(ceil(length(spiketrain)/15),1);
    for n = 1:ceil(length(spiketrain)/15)
        if n ~= ceil(length(spiketrain)/15)
            rebinned_spiketrain(n) = sum(spiketrain(((n-1)*15)+1:((n-1)*15)+15));
        else
            rebinned_spiketrain(n) = sum(spiketrain(((n-1)*15)+1:end));
        end
    end
    
    velx = diff([posx_c(1); posx_c]); vely = diff([posy_c(1); posy_c]);
    speed = sqrt(velx.^2+vely.^2)*sampleRate; 
    maxSpeed = 50; speed(speed>maxSpeed) = maxSpeed;
    
    speed_scores(p) = corr(speed,spiketrain); % speed score of cell
    
    for i = 1:100
        min_shift = ceil(20/dt); % min shift is 20 s
        max_shift = length(spiketrain)-(20/dt); % max shift is length of trial minus 20 s 
        rand_shift = round(min_shift + rand*(max_shift-min_shift)); % amount to shift spiketrain by
        
        spiketrain_shift = circshift(spiketrain,rand_shift); % shifted spiketrain
        
        idx = ((p-1)*100)+i;
        shifted_speed_scores(idx,1) = p; % cell num
        shifted_speed_scores(idx,2) = corr(speed,spiketrain_shift); % Pearson correlation between speed and shifted spiketrain
    end
end

sig_percentile = prctile(shifted_speed_scores(:,2),99);
% speed_cells_idx = find(shifted_speed_scores(:,2) > sig_percentile);
% speed_cells = unique(shifted_speed_scores(speed_cells_idx,1));
speed_cells = find(speed_scores > sig_percentile);

save speed_cell_data.mat
