function [fourierSpectrogram, polarSpectrogram, smoothRho, superimpose, beforeMaxPower, maxPower, nComponents, isPeriodic] = fourier_transform(rateMap, meanFr, spiketrain, dt, posx, posy)
% Calculates the Fourier Spectrogram of a given rate map and determines the
% main components and polar distribution

% Adapted from Neural Representations of Location Composed of Spatially 
% Periodic Bands, Krupic et al 2012

% Performs 2D FFT to calculate the fourier spectrogram and max power for
% the cell
fourierSpectrogram = fft2(rateMap,256,256);
beforeMaxPower = max(abs(fftshift(fourierSpectrogram(:))));

fourierSpectrogram = fourierSpectrogram ./ (meanFr * sqrt(size(rateMap, 1) * size(rateMap, 2)));
fourierSpectrogram = abs(fftshift(fourierSpectrogram)); % Matrix values represent power of Fourier spectrum

maxPower = max(fourierSpectrogram(:));

% Shuffles the spiketrain of the cell 150 times and calculates the max
% power of the shuffled fourier spectrogram. 
shiftedMaxPowers = zeros(150,1);
for i = 1:150
    minShift = ceil(20/dt); % min shift is 20 s
    maxShift = length(spiketrain)-(20/dt); % max shift is length of trial minus 20 s 
    randShift = round(minShift + rand*(maxShift-minShift)); % amount to shift spiketrain by
        
    shiftedSpiketrain = circshift(spiketrain,randShift); % shifted spiketrain
    shiftedFiringRate = shiftedSpiketrain ./ dt; 
        
    shiftedRateMap = compute_2d_tuning_curve(posx,posy,shiftedFiringRate,20,0,100);
    shiftedRateMap = shiftedRateMap - nanmean(shiftedRateMap(:));
    
    shiftedRateMap(isnan(shiftedRateMap)) = 0; % Converts all NaNs to zeros for FFT calculation
    
    shiftedFourier = fft2(shiftedRateMap,256,256);
    shiftedFourier = shiftedFourier ./ (meanFr * sqrt(size(shiftedRateMap, 1) * size(shiftedRateMap, 2)));
    shiftedFourier = abs(fftshift(shiftedFourier)); % Matrix values represent power of Fourier spectrum
    shiftedMaxPowers(i) = max(shiftedFourier(:));
end

% Determines if the cell is periodic. A cell is periodic if its max power
% is greater than the 95th percentile of shuffled data. 
sigPercentile = prctile(shiftedMaxPowers, 95);
if (maxPower > sigPercentile)
    isPeriodic = true; 
else
    isPeriodic = false;
end

% Reduces effects of noise by subtracting the 50th percentile value of
% power from the Fourier spectrogram
threshold1 = prctile(shiftedMaxPowers, 50);
polarSpectrogram = fourierSpectrogram - threshold1;
polarSpectrogram(polarSpectrogram < 0) = 0;

% Excludes spatial sub-harmonic frequencies
threshold2 = max(polarSpectrogram(:)) * 0.1;
mainComponents = polarSpectrogram;
mainComponents(mainComponents <= threshold2) = 0;

% Calculates the average power along every orientation of the polar spectrogram 
rhoMatrix = zeros(180, 2);
% rhoMatrix = zeros(361,2);
% only look at top half because of symmetry
for j = 1:length(mainComponents) / 2 % row
    for k = 1:length(mainComponents) % col
        y = 128.5 - j;
        x = k - 128.5;
        if x > 0
            if y > 0
                theta = floor(rad2deg(atan(y/x)));
            else
                theta = floor(rad2deg(atan(y/x) + 2*pi));
            end
        else
            theta = floor(rad2deg(atan(y/x) + pi));
        end
        rho = mainComponents(j,k);
        if (rho ~= 0 && ~isnan(rho))
            rhoMatrix(theta+1,1) = rhoMatrix(theta+1,1) + rho;
            rhoMatrix(theta+1,2) = rhoMatrix(theta+1,2) + 1;
        end
    end
end

rhoMeanPower = rhoMatrix(:,1) ./ rhoMatrix(:,2);
rhoMeanPower(isnan(rhoMeanPower)) = 0;

% Smooth curve
gaussFilter = gausswin(17);
gaussFilter = gaussFilter / sum(gaussFilter);
smoothRho = conv(rhoMeanPower, gaussFilter,'same');
smoothRho = [smoothRho; smoothRho; smoothRho(1)]; 

% Finds peaks and number of components of a polar Fourier spectrogram. 
% Doesn't count peaks within 40 degrees of a larger peak. 
[peaks, locs] = findpeaks(smoothRho);
mainPeaks = []; mainLocs = [];
for loc = locs'
    if (loc > 180)
        continue
    end
    
    % Determines range of angles to look for neighboring peaks in
    % Wraps angles between 0 and 360
    minAngle = loc - 20;
    if minAngle < 0
        minAngle = minAngle + 360;
    end
    maxAngle = loc + 20; 
    if maxAngle > 360
        maxAngle = maxAngle - 360;
    end
    
    if (minAngle < maxAngle)
        neighborPeaks = peaks(locs >= minAngle & locs <= maxAngle);
    else
        neighborPeaks = peaks((locs >= minAngle & locs <= 360) | (locs >= 0 & locs <= maxAngle));
    end
    
    % Only count highest peak in 40 degree range
    highestPeak = max(neighborPeaks);
    highestLoc = max(locs(peaks == highestPeak));
    
    if (highestLoc > 180)
        highestLoc = highestLoc - 180;
    end
    
    mainPeaks = [mainPeaks highestPeak];
    mainLocs = [mainLocs highestLoc];
end
mainPeaks = unique(mainPeaks);
mainLocs = unique(mainLocs);
nComponents = length(mainLocs);
mainPeaks(mainPeaks ~= 0) = ceil(max(mainPeaks));

% Stores polar graph data of the polar Fourier Spectrogram's power
superimposeTheta = kron(mainLocs - 1, [0 1]);
superimposeRho = kron(mainPeaks, [0 1]);
superimpose = cell(2,1);
superimpose{1} = superimposeTheta;
superimpose{2} = superimposeRho;

end

