function [maxCircularField, maxGridScore] = find_circular_field(rateMap)
% Determines the circular region surrounding multiple firing fields to 
% calculate the grid score from given a rate map 

% git test

% Threshold firing rate for a position bin to be considered a peak
maxRate = max(rateMap(:));
threshold = 0.2 * maxRate;

% Converts all position bins with firing rates below the threshold to zero
% and labels connected components with firing rates above the threshold
modifiedRateMap = rateMap;
modifiedRateMap(modifiedRateMap <= threshold | isnan(modifiedRateMap)) = 0;
modifiedRateMap(modifiedRateMap > threshold) = 1;
modifiedRateMap = bwlabel(modifiedRateMap);
        
% Increases number of pixels per time bin to enhance image and circle quality/accuracy
expand = 3;
modifiedRateMapExpanded = kron(modifiedRateMap, ones(expand)); % Map of peaks 
rateMapExpanded = kron(rateMap, ones(expand)); % Original map of all firing rates

% Loops through all peaks and centers a circular field around that peak.
% Chooses the circular field with the highest grid score. 
maxGridScore = -inf;
maxCircularField = zeros(size(rateMapExpanded));
for i = 1:max(modifiedRateMap(:))
    % Only considers connected components with greater than 6*expand
    % position bins as peaks 
    nBins = length(find(modifiedRateMapExpanded == i));
    if (nBins > 6 * expand)
        % Gets coordinates of the center peak to calculate the center of
        % the circle
        [row,col] = find(modifiedRateMapExpanded == i);
        xCenter = ceil((max(col) + min(col)) / 2);
        yCenter = ceil((max(row) + min(row)) / 2);
        center = [xCenter yCenter];
            
        % Gets coordinates of all peaks besides the center peak
        [row1, col1] = find(modifiedRateMapExpanded ~= i & modifiedRateMapExpanded ~= 0);
        coordinates = [col1 row1];    
        distances = sqrt(sum(bsxfun(@minus, coordinates, center).^2,2));
        % Number of occurances of each element in the modified rate map
        temp = modifiedRateMapExpanded;
        temp(temp == 0) = [];
        histogram = hist(temp(:));
        maxBins = max(histogram);
        % Radius of circle equals the minimum distance between the
        % center and another peak plus the approxmate diameter of the
        % largest peak in the rate map 
        radius = ceil(min(distances)) + ceil(sqrt(maxBins));
        
        % Extracts the circular field from the rate map by multiplying a
        % mask of 0s and 1s with the rate map 
        dim = size(rateMapExpanded);
        yDim = dim(1); xDim = dim(2);
        [xMask,yMask] = meshgrid(-(xCenter-1):(xDim-xCenter),-(yCenter-1):(yDim-yCenter));
        mask = ((xMask.^2 + yMask.^2) <= radius^2);
        circularField = rateMapExpanded .* mask;
           
        % Concatanates zero vectors horizontally and vertically to center the
        % circular field for later crosscorrelation calculations
        horizontalShift = abs(xCenter - xDim + xCenter);
        verticalShift = abs(yCenter - yDim + yCenter);
        dimensions = size(rateMapExpanded);
        height = dimensions(1);
        width = dimensions(2);
        horizontalAdd = zeros(height, horizontalShift);
        verticalAdd = zeros(verticalShift, width + horizontalShift);
        if (xCenter > length(rateMapExpanded) / 2)
            circularField = horzcat(circularField, horizontalAdd);
        else
            circularField = horzcat(horizontalAdd, circularField);
        end
            
        if (yCenter > length(rateMapExpanded) / 2)
            circularField = vertcat(circularField, verticalAdd);
        else
            circularField = vertcat(verticalAdd, circularField);
        end
            
        gridScore = calculate_grid_score(circularField);
        if (gridScore > maxGridScore)
            maxGridScore = gridScore;
            maxCircularField = circularField;
        end
    end
end

end
