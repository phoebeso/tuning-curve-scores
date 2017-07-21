clear all; clc; 

% load spikePos_openfield_D4.mat
% load thetaEEG_JZ1_day4_nt11.mat
% load JZ1posday4.mat
% posStruct6 = pos{1,4}{1,6};
% posStruct8 = pos{1,4}{1,8};
% spikePosStruct = JZ1spikesPosD4;
% thetaStruct6 = thetaEEG_JZ1_day4_epoch6_nt11;
% thetaStruct8 = thetaEEG_JZ1_day4_epoch8_nt11;
% 
% nCells = length(spikePosStruct); % number of cells being analyzed
% post = posStruct6.data(:,1);
% posx = posStruct6.data(:,2);
% posy = posStruct6.data(:,3);
% filt_eeg = thetaStruct6.data;
% 
% posx = posx - min(posx);
% posy = posy - min(posy);
% 
% nPosBins = 61; % each pos bin is 2x2 cm 
% boxSize = ceil(max(max(posx),max(posy)));
% 
% gridScores = zeros(nCells,1);
% 
% for n = 52:52
% 
%     spiketimes = spikePosStruct(n).spikes(:,1);
%     timebins = linspace(post(1),post(end)+post(3)-post(2),length(post)+1);
%     spiketrain = histcounts(spiketimes,timebins)';
% 
%     filter = gaussian(-4:4,2,0); filter = filter / sum(filter); 
%     dt = post(3) - post(2); firingRate = spiketrain / dt;
%     smoothFiringRate = conv(firingRate,filter,'same');
%    
%     
%     [rateMap] = compute_2d_tuning_curve(posx,posy,smoothFiringRate,nPosBins,0,boxSize);
%     gridScores(n) = calculate_grid_score(firingFields);

%     [autocorrelogram] = compute_autocorrelogram(firingFields);
    
%     figure(1)
%     imagesc(rateMap); colorbar
    
%     figure(2)
%     imagesc(firingFields, [floor(min(rateMap(:))),ceil(max(rateMap(:)))]); colorbar
%     figure(3)
%     imagesc(autocorrelogram); colorbar
% end

% sigPercentile = prctile(gridScores,99);
% gridCells = find(gridScores > sigPercentile);

load data_for_cell77.mat

posx = posx_c;
posy = posy_c;
nPosBins = 20;
boxSize = 100;

filter = gaussian(-4:4,2,0); filter = filter / sum(filter); 
dt = post(3) - post(2); firingRate = spiketrain / dt;
smoothFiringRate = conv(firingRate,filter,'same');

[rateMap] = compute_2d_tuning_curve(posx,posy,smoothFiringRate,nPosBins,0,boxSize);
[circularField, gridScore] = find_circular_field(rateMap);

figure(1)
imagesc(rateMap); colorbar
figure(2); colorbar
imagesc(circularField, [floor(min(rateMap(:))),ceil(max(rateMap(:)))]); colorbar
