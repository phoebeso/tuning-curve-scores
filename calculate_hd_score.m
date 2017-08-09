function [hdOccupancy, hdRates, hdScore] = calculate_hd_score(direction,spiketrain,dt,nHdBins)
% Calculates the head direction score and tuning curve. Head direction
% tuning curve calculated by first binning the head directions into 6
% degree bins. 

% calculating, for the entire trial, the length of the mean vector of firing 
% rate as a function of direction (directional tuning).
% http://www.sciencedirect.com/science/article/pii/S0960982213015212

hdBins = (0:2*pi/nHdBins:2*pi-2*pi/nHdBins)'; % Represent left edge of bin
hdCount = zeros(nHdBins,2); % col 1 = num of times hd bin accessed, col 2 = num of spikes total

for i = 1:numel(direction)
   
    % figure out the hd index
    diference = direction(i) - hdBins;
    [~, idx] = min(diference(diference >= 0));
    hdCount(idx,1) = hdCount(idx,1) + 1;
    hdCount(idx,2) = hdCount(idx,2) + spiketrain(i);

end

hdOccupancy = hdCount(:,1);
hdRates = hdCount(:,2) ./ (hdCount(:,1) .* dt); % vector with spikerates for each hd bin

% convert coordinates from polar to rectangular form
x = hdRates .* cos(hdBins);
y = hdRates .* sin(hdBins);
sumLength = sum(sqrt(x.^2 + y.^2));

% computes sum vector
rx = sum(x);
ry = sum(y);
rLength = sqrt(rx^2 + ry^2);

hdScore = rLength / sumLength; 

end
