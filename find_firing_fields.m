function [maxFiringFields, maxGridScore] = find_firing_fields(rateMap)
% Determines the circular region surrounding multiple firing fields given a
% rate map 

peakRate = max(rateMap(:));
threshold = 0.2 * peakRate;

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

nPeaks = numel(unique(modifiedRateMap)) - 1;
        
% Increases number of pixels per time bin to enhance image and circle quality/accuracy 
modifiedRateMapExpanded = kron(modifiedRateMap, ones(3));
rateMapExpanded = kron(rateMap, ones(3));

maxGridScore = -inf;
maxFiringFields = zeros(size(rateMapExpanded));
for i = 1:nPeaks
    [row,col] = find(modifiedRateMapExpanded == i);
    xCenter = ceil((max(col) + min(col)) / 2);
    yCenter = ceil((max(row) + min(row)) / 2);
    center = [xCenter yCenter];
    
    [row1, col1] = find(modifiedRateMapExpanded ~= i & modifiedRateMapExpanded ~= 0);
    points = [col1 row1];
    
    distances = sqrt(sum(bsxfun(@minus, points, center).^2,2));
    histogram = sort(hist(modifiedRateMapExpanded(:)), 'descend'); % number of occurances of each element in the modified rate map
    n = histogram(2); % gets second largest value, because zero (no peak) should in theory be the value with the highest number of occurances
    
    radius = ceil(min(distances)) + ceil(sqrt(n)); % radius is the min distance from the center to another peak plus the approximate diameter of the
    % largest peak in the rate map 
    
    dim = size(modifiedRateMapExpanded);
    xDim = dim(1); yDim = dim(2);
    [xx,yy] = ndgrid((1:yDim)-yCenter,(1:xDim)-xCenter);
    mask = zeros(xDim,yDim);
    mask((xx.^2 + yy.^2) < radius^2) = 1;
    firingFields = rateMapExpanded .* mask;
    
%     % Center the circle for later crosscorrelation calculations
%     xShift = ceil(xDim/2) - xCenter;
%     yShift = ceil(yDim/2) - yCenter;
%     firingFields = circshift(firingFields,[yShift,xShift]);

    % Concatanates zero vectors horizontally and vertically to center the
    % firing field for later crosscorrelation calculations 
    horizontalShift = abs(xCenter - length(rateMapExpanded) + xCenter);
    verticalShift = abs(yCenter - length(rateMapExpanded) + yCenter);
    dimensions = size(rateMapExpanded);
    height = dimensions(1);
    width = dimensions(2);
    horizontalAdd = zeros(height, horizontalShift);
    verticalAdd = zeros(verticalShift, width + horizontalShift);
    if (xCenter > length(rateMapExpanded) / 2)
        firingFields = horzcat(firingFields, horizontalAdd);
    else
        firingFields = horzcat(horizontalAdd, firingFields);
    end
    
    if (yCenter > length(rateMapExpanded) / 2)
        firingFields = vertcat(firingFields, verticalAdd);
    else
        firingFields = vertcat(verticalAdd, firingFields);
    end
    
    gridScore = calculate_grid_score(firingFields);
    if (gridScore > maxGridScore)
        maxGridScore = gridScore;
        maxFiringFields = firingFields;
    end
end

% Gets coordinates of firing fields 
% [row,col] = find(modifiedRateMapExpanded ~= 0);
% coordinates = [col row];
%     
% % Fits circle around peak coordinates
% circleParams = circle_fit_by_pratt(coordinates);
% xCenter = ceil(circleParams(1));
% yCenter = ceil(circleParams(2));
% radius = ceil(circleParams(3)) + 3;

% Gets circular area fitted around peaks from rate map 
% dim = size(modifiedRateMapExpanded);
% xDim = dim(1); yDim = dim(2);
% [xx,yy] = ndgrid((1:yDim)-yCenter,(1:xDim)-xCenter);
% mask = zeros(xDim,yDim);
% mask((xx.^2 + yy.^2) < radius^2) = 1;
% firingFields = rateMapExpanded .* mask;
%     
% % Center the circle for later crosscorrelation calculations 
% xShift = ceil(xDim/2) - xCenter;
% yShift = ceil(yDim/2) - yCenter;
% firingFields = circshift(firingFields,[yShift,xShift]);
%     
% figure(1)
% imagesc(firingFields); colorbar 

end
