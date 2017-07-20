clear all; clc; 

load spikePos_openfield_D4.mat
load thetaEEG_JZ1_day4_nt11.mat
load JZ1posday4.mat
posStruct6 = pos{1,4}{1,6};
posStruct8 = pos{1,4}{1,8};
spikePosStruct = JZ1spikesPosD4;
thetaStruct6 = thetaEEG_JZ1_day4_epoch6_nt11;
thetaStruct8 = thetaEEG_JZ1_day4_epoch8_nt11;
nCells = length(spikePosStruct); % number of cells being analyzed

cellNum = 32;
post = posStruct6.data(:,1);
posx = posStruct6.data(:,2);
posy = posStruct6.data(:,3);
filt_eeg = thetaStruct6.data;

posx = posx - min(posx);
posy = posy - min(posy);

nPosBins = 61; % each pos bin is 2x2 cm 
boxSize = ceil(max(max(posx),max(posy)));

spiketimes = spikePosStruct(cellNum).spikes(:,1);
timebins = linspace(post(1),post(end)+post(3)-post(2),length(post)+1);
spiketrain = histcounts(spiketimes,timebins)';

filter = gaussian(-4:4,2,0); filter = filter / sum(filter); 
dt = post(3) - post(2); firingRate = spiketrain / dt;
smoothFiringRate = conv(firingRate,filter,'same');

[rateMap] = compute_2d_tuning_curve(posx,posy,smoothFiringRate,nPosBins,0,boxSize);

% original rate map
% pos_curve_expand = kron(rate_map, ones(3));
% figure(1);
% imagesc(pos_curve_expand); colorbar
% rotated rate map
% figure(2)
% angle = 60;
% rotated_img = imrotate(pos_curve_expand, angle, 'crop'); 
% rotated_img(rotated_img == 0) = NaN;
% imagesc(rotated_img); colorbar

[firingFields] = find_firing_fields(rateMap);
[gridScore] = calculate_grid_score(firingFields);

% make mask and crop circular area from image
% x_center = 30;
% y_center = 50;
% radius = 40;
% dim = size(pos_curve_expand);
% x_dim = dim(1);
% y_dim = dim(2);
% [xx,yy] = ndgrid((1:y_dim)-y_center,(1:x_dim)-x_center); % 125 = dimensions, 50, 100 = y-coordinate, x-coordinate
% mask = zeros(x_dim, y_dim);
% mask((xx.^2 + yy.^2) < radius^2) = 1; % 25 = radius of circle
% circle_graph = pos_curve_expand.*mask;
% % center circle
% x_shift = ceil(x_dim/2) - x_center;
% y_shift = ceil(y_dim/2) - y_center;
% circle_graph = circshift(circle_graph,[y_shift,x_shift]);
% figure(3)
% imagesc(circle_graph, [floor(min(rate_map(:))),ceil(max(rate_map(:)))]); colorbar
% % rotate circle
% figure(4)
% rotated_circle_graph = imrotate(circle_graph,angle,'crop');
% rotated_circle_graph(rotated_circle_graph == 0) = NaN;
% imagesc(rotated_circle_graph, [floor(min(rate_map(:))),ceil(max(rate_map(:)))]); colorbar
% 
% crosscorrelation = calculate_crosscorrelation(circle_graph, rotated_circle_graph);
