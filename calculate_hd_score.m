function [hdRates, hdScore] = calculate_hd_score(direction,spiketrain,dt,nHdBins)
% Calculates the head direction score and tuning curve. Head direction
% tuning curve calculated by first binning the head directions into 6
% degree bins. 

hdBins = 0:2*pi/nHdBins:2*pi-2*pi/nHdBins;
hdCount = zeros(nHdBins,2); % col 1 = num of times hd bin accessed, col 2 = num of spikes total

for i = 1:numel(direction)
   
    % figure out the hd index
    [~, idx] = min(abs(direction(i) - hdBins));
    hdCount(idx,1) = hdCount(idx,1) + 1;
    hdCount(idx,2) = hdCount(idx,2) + spiketrain(i);

end
    
hdRates = hdCount(:,2) ./ (hdCount(:,1).*dt); % vector with spikerates for each hd bin 
% smooth = hamming(10);
% smoothHdRates = conv(hdRates,smooth,'same'); % Smooth head direction tuning curve for display purposes

hdScore = circ_r(hdRates, [], 2*pi/nHdBins, []);

end
