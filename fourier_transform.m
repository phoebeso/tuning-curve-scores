function [fourierSpectrogram, polarSpectrogram, nComponents] = fourier_transform(matrix, meanFr, spiketrain, dt, posx, posy)

fourierSpectrogram = fft2(matrix,256,256);
fourierSpectrogram = abs(fftshift(fourierSpectrogram)); % Matrix values represent power of Fourier spectrum
fourierSpectrogram = fourierSpectrogram ./ (meanFr * length(matrix));

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
    fourierSpectrogram2 = fft2(shiftedRateMap,256,256);
    fourierSpectrogram2 = abs(fftshift(fourierSpectrogram2)); % Matrix values represent power of Fourier spectrum
    fourierSpectrogram2 = fourierSpectrogram2 ./ (meanFr * length(shiftedRateMap));
    shiftedPower(i) = max(fourierSpectrogram2(:));
end

threshold1 = prctile(shiftedPower, 50);

% Reduces effects of noise by subtracting the 50th percentile value of
% power from the Fourier spectrogram
% reshapeMatrix = reshape(fourierSpectrogram, [], 1);
% % threshold1 = prctile(reshapeMatrix, 50);
% threshold1 = maxPower * 0.5;
polarSpectrogram = fourierSpectrogram - threshold1;
polarSpectrogram(polarSpectrogram < 0) = 0;

% Excludes spatial sub-harmonic frequencies
threshold2 = max(polarSpectrogram(:)) * 0.1;
modifiedMatrix = polarSpectrogram;
modifiedMatrix(modifiedMatrix <= threshold2) = 0;
modifiedMatrix(modifiedMatrix > threshold2) = 1;
modifiedMatrix = bwlabel(modifiedMatrix); % IDs the peaks 

nComponents = max(modifiedMatrix(:));

end

