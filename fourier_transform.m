function [fourierSpectrogram, polarSpectrogram, smoothRho, superimpose, maxPower, nComponents] = fourier_transform(rateMap, meanFr, spiketrain, dt, posx, posy)
% Calculates the Fourier Spectrogram of a given rate map and determines the
% main components and polar distribution

% Adapted from Neural Representations of Location Composed of Spatially 
% Periodic Bands, Krupic et al 2012

fourierSpectrogram = fft2(rateMap,256,256);
fourierSpectrogram = fourierSpectrogram ./ (meanFr * sqrt(size(rateMap, 1) * size(rateMap, 2)));
fourierSpectrogram = abs(fftshift(fourierSpectrogram)); % Matrix values represent power of Fourier spectrum

maxPower = max(fourierSpectrogram(:));

shiftedPower = zeros(100,1);
for i = 1:100
    minShift = ceil(20/dt); % min shift is 20 s
    maxShift = length(spiketrain)-(20/dt); % max shift is length of trial minus 20 s 
    randShift = round(minShift + rand*(maxShift-minShift)); % amount to shift spiketrain by
        
    shiftedSpiketrain = circshift(spiketrain,randShift); % shifted spiketrain
    shiftedFiringRate = shiftedSpiketrain ./ dt; 
        
    shiftedRateMap = compute_2d_tuning_curve(posx,posy,shiftedFiringRate,20,0,100);
    shiftedRateMap(isnan(shiftedRateMap)) = 0; 
    shiftedRateMap = shiftedRateMap - mean(shiftedRateMap(:));
    shiftedFourier = fft2(shiftedRateMap,256,256);
    shiftedFourier = shiftedFourier ./ (meanFr * sqrt(size(shiftedRateMap, 1) * size(shiftedRateMap, 2)));
    shiftedFourier = abs(fftshift(shiftedFourier)); % Matrix values represent power of Fourier spectrum
    shiftedPower(i) = max(shiftedFourier(:));
end

% Reduces effects of noise by subtracting the 50th percentile value of
% power from the Fourier spectrogram
threshold1 = prctile(shiftedPower, 50);
polarSpectrogram = fourierSpectrogram - threshold1;
polarSpectrogram(polarSpectrogram < 0) = 0;

% Excludes spatial sub-harmonic frequencies
threshold2 = max(polarSpectrogram(:)) * 0.1;
mainComponents = polarSpectrogram;
mainComponents(mainComponents <= threshold2) = 0;

% Calculates the average power along every orientation of the polar spectrogram 
rhoMatrix = zeros(361, 2);
for j = 1:length(mainComponents) % row
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
        if rho ~= 0
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

% Finds peaks in polar Fourier spectrogram and don't count peaks within 10 
% degrees of a larger peak 
[peaks, locs] = findpeaks(smoothRho);
mainPeaks = []; mainLocs = [];
for loc = locs'
    % Determines range of angles to look for neighboring peaks in
    % Wraps angles between 0 and 360
    minAngle = loc - 10;
    if minAngle < 0
        minAngle = minAngle + 360;
    end
    maxAngle = loc + 10; 
    if maxAngle > 360
        maxAngle = maxAngle - 360;
    end
    
    if (minAngle < maxAngle)
        neighborPeaks = peaks(locs >= minAngle & locs <= maxAngle);
    else
        neighborPeaks = peaks((locs >= minAngle & locs <= 360) | (locs >= 0 & locs <= maxAngle));
    end
    
    % Only count highest peak in 10 degree range
    highestPeak = max(neighborPeaks);
    highestLoc = locs(peaks == highestPeak);
    
    if (highestLoc > 180)
        continue
    end
    
    mainPeaks = [mainPeaks highestPeak];
    mainLocs = [mainLocs highestLoc];
end
nComponents = length(mainPeaks);

superimposeTheta = kron(mainLocs - 1, [0 1]);
superimposeRho = kron(mainPeaks, [0 1]);
superimpose = cell(2,1);
superimpose{1} = superimposeTheta;
superimpose{2} = superimposeRho;

end

