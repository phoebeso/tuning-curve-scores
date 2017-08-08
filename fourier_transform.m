function [shiftedSpectrogram, polarSpectrogram, nComponents] = fourier_transform(correlationMatrix)

fourierSpectrogram = fft2(correlationMatrix,256,256);
shiftedSpectrogram = abs(fftshift(fourierSpectrogram));

maxPower = max(shiftedSpectrogram(:));
polarSpectrogram = shiftedSpectrogram - (maxPower * 0.5);
polarSpectrogram(polarSpectrogram < 0) = 0;

threshold = max(polarSpectrogram(:)) * 0.1;
modifiedMatrix = polarSpectrogram;
modifiedMatrix(modifiedMatrix <= threshold | isnan(modifiedMatrix)) = 0;
modifiedMatrix(modifiedMatrix > threshold) = 1;
modifiedMatrix = bwlabel(modifiedMatrix); % IDs the peaks 

nComponents = max(modifiedMatrix(:));

end

