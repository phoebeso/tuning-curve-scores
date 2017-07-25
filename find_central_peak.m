function [maxField, maxGridScore] = find_central_peak(rateMap)
% Determines the circular region surrounding multiple firing fields to 
% calculate the grid score from given rate map 
 
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
maxField = zeros(size(rateMapExpanded));
for i = 1:max(modifiedRateMap(:))
    % Only considers connected components with greater than 6*expand
    % position bins as peaks 
    peakSize = length(find(modifiedRateMapExpanded == i));
    if (peakSize > 6 * expand)
        % Gets coordinates of the center peak to calculate the center of
        % rotation
        [rowCenter, colCenter] = find(modifiedRateMapExpanded == i);
        xCenter = ceil((max(colCenter) + min(colCenter)) / 2);
        yCenter = ceil((max(rowCenter) + min(rowCenter)) / 2);
           
        % Concatanates nan vectors horizontally and vertically to center the
        % circular field for later crosscorrelation calculations
        % Additionally shifts the central peak's rows and cols for
        % identification purposes for calculations
        dim = size(rateMapExpanded);
        yDim = dim(1); xDim = dim(2);
        horizontalShift = abs(xCenter - xDim + xCenter);
        verticalShift = abs(yCenter - yDim + yCenter);
        horizontalAdd = nan(yDim, horizontalShift);
        verticalAdd = nan(verticalShift, xDim + horizontalShift);
        shiftedRateMap = rateMapExpanded; 
        if (xCenter > length(rateMapExpanded) / 2)
            shiftedRateMap = horzcat(shiftedRateMap, horizontalAdd);
            shiftedColCenter = colCenter;
        else
            shiftedRateMap = horzcat(horizontalAdd, shiftedRateMap);
            shiftedColCenter = colCenter + horizontalShift; 
        end
            
        if (yCenter > length(rateMapExpanded) / 2)
            shiftedRateMap = vertcat(shiftedRateMap, verticalAdd);
            shiftedRowCenter = rowCenter;
        else
            shiftedRateMap = vertcat(verticalAdd, shiftedRateMap);
            shiftedRowCenter = rowCenter + verticalShift; 
        end
        
        % Excludes the central peak from the grid score calculation 
        fieldWithoutCenter = shiftedRateMap;
        fieldWithoutCenter(shiftedRowCenter, shiftedColCenter) = 0;
        gridScore = calculate_grid_score(fieldWithoutCenter);
        if (gridScore > maxGridScore)
            maxGridScore = gridScore;
            maxField = shiftedRateMap;
        end
    end 
end

end
