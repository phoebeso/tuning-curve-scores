function [firingFields] = find_firing_fields(rateMap)
% Determines the circular region surrounding multiple firing fields given a
% rate map 

peakRate = max(rateMap(:));
threshold = 0.8 * peakRate;

% Converts all position bins with firing rates below the threshold to zero
% and labels connected components with firing rates above the threshold
modifiedRateMap = rateMap;
modifiedRateMap(modifiedRateMap <= threshold | isnan(modifiedRateMap)) = 0;
modifiedRateMap(modifiedRateMap > threshold) = 1;
modifiedRateMap = bwlabel(modifiedRateMap);
    
% Removes connected componenets with less than 6 position bins 
for i = 1:max(modifiedRateMap(:))
    nPosBins = find(modifiedRateMap == i);
    if (nPosBins < 6)
        modifiedRateMap(modifiedRateMap == i) = 0;
    end
end
    
%     [radius, center] = ExactMinBoundCircle(coordinates); % radius and center of original circle 
    
% Increases number of pixels per time bin to enhance image and circle quality 
modifiedRateMapExpanded = kron(modifiedRateMap, ones(3));
rateMapExpanded = kron(rateMap, ones(3));

% Gets coordinates of firing fields 
[row,col] = find(modifiedRateMapExpanded ~= 0);
coordinates = [col row];
    
%     [radius, center] = ExactMinBoundCircle(coordinates);
%     [x_center, y_center, radius] = circfit(col,row);
%     x_center = ceil(x_center);
%     y_center = ceil(y_center);

% Fits circle around peak coordinates
circleParams = CircleFitByPratt(coordinates);
xCenter = ceil(circleParams(1));
yCenter = ceil(circleParams(2));
radius = ceil(circleParams(3)) + 3;

% Gets circular firing field from rate map fitted around peaks 
dim = size(modifiedRateMapExpanded);
%     x_center = (center(1) * 3) - 1;
%     y_center = (center(2) * 3) - 1;
%     radius = (ceil(radius) * 3) + 3;
%     x_center = round(center(1));
%     y_center = round(center(2));
xDim = dim(1); yDim = dim(2);
[xx,yy] = ndgrid((1:yDim)-yCenter,(1:xDim)-xCenter);
mask = zeros(xDim, yDim);
mask((xx.^2 + yy.^2) < radius^2) = 1;
firingFields = rateMapExpanded .* mask;
    
% Center the circle for later crosscorrelation calculations 
x_shift = ceil(x_dim/2) - x_center;
y_shift = ceil(y_dim/2) - y_center;
firingFields = circshift(firingFields,[y_shift,x_shift]);
    
figure(1)
imagesc(firingFields); colorbar 

end
