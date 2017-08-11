function [fourierSpectrogram, polarSpectrogram, nComponents] = fourier_transform(rateMap, meanFr, spiketrain, dt, posx, posy)

fourierSpectrogram = fft2(rateMap,256,256);
fourierSpectrogram = abs(fftshift(fourierSpectrogram)); % Matrix values represent power of Fourier spectrum
fourierSpectrogram = fourierSpectrogram ./ (meanFr * sqrt(size(rateMap, 1) * size(rateMap, 2)));

maxPower = max(fourierSpectrogram(:));

shiftedPower = zeros(100,1);
for i = 1:100
    minShift = ceil(20/dt); % min shift is 20 s
    maxShift = length(spiketrain)-(20/dt); % max shift is length of trial minus 20 s 
    randShift = round(minShift + rand*(maxShift-minShift)); % amount to shift spiketrain by
        
    shiftedSpiketrain = circshift(spiketrain,randShift); % shifted spiketrain
    shiftedFiringRate = shiftedSpiketrain ./ dt; 
        
    shiftedRateMap = compute_2d_tuning_curve(posx,posy,shiftedFiringRate,20,0,100);
    shiftedRateMap = shiftedRateMap - mean(shiftedRateMap(:));
    shiftedFourier = fft2(shiftedRateMap,256,256);
    shiftedFourier = abs(fftshift(shiftedFourier)); % Matrix values represent power of Fourier spectrum
    shiftedFourier = shiftedFourier ./ (meanFr * sqrt(size(shiftedRateMap, 1) * size(shiftedRateMap, 2)));
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

% mainComponents(mainComponents > threshold2) = 1;
% mainComponents = bwlabel(mainComponents); % IDs the peaks 

% nComponents = max(mainComponents(:));

rhoArray = zeros(361, 1);
for j = 1:length(mainComponents) % row
    for k = 1:length(mainComponents) % col
        y = 128.5 - j;
        x = k - 128.5;
        if x > 0
            if y > 0
                theta = round(rad2deg(atan(y/x)));
            else
                theta = round(rad2deg(atan(y/x) + 2*pi));
            end
        else
            theta = round(rad2deg(atan(y/x) + pi));
        end
        rho = mainComponents(j,k);
        rhoArray(theta + 1) = max(rhoArray(theta+1), rho);
    end
end

angles = 0:360;

end

