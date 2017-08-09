function [hdOccupancy, hdRates, hdScore, xmean, ymean] = calculate_hd_score(direction,spiketrain,dt,nHdBins)
% Calculates the head direction score and tuning curve. Head direction
% tuning curve calculated by first binning the head directions into 6
% degree bins. 

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

normalizedHdRates = hdRates ./ max(hdRates); 

% convert coordinates from polar to rectangular form
x = normalizedHdRates .* cos(hdBins);
y = normalizedHdRates .* sin(hdBins);

% obtain length
xmean = sum((normalizedHdRates ./ max(normalizedHdRates)) .* x) / sum(x ~= 0);
ymean = sum((normalizedHdRates ./ max(normalizedHdRates)) .* y) / sum(y ~= 0);
xmean = sum(x) / sum(x ~= 0);
ymean = sum(y) / sum(y ~= 0);
xmean = sum(x);
ymean = sum(y);
hdScore = sqrt(xmean^2 + ymean^2);

end
