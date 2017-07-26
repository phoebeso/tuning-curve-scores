function [maxShiftedField, maxGridScore] = find_central_peak(rateMap)
% Given a rate map, determines the center peak to rotate the map around
% from multiple firing fields that maximizes the grid score 

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
maxShiftedField = rateMapExpanded;
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
           
        % Concatanates zero vectors horizontally and vertically to center the
        % circular field for later crosscorrelation calculations
        % Additionally shifts the central peak's rows and cols for
        % identification purposes for calculations
        dim = size(rateMapExpanded);
        yDim = dim(1); xDim = dim(2);
        horizontalShift = abs(xCenter - xDim + xCenter);
        verticalShift = abs(yCenter - yDim + yCenter);
        horizontalAdd = zeros(yDim, horizontalShift);
        verticalAdd = zeros(verticalShift, xDim + horizontalShift);
        if (xCenter > length(rateMapExpanded) / 2)
            shiftedRateMap = [rateMapExpanded horizontalAdd];
            shiftedColCenter = colCenter;
        else
            shiftedRateMap = [horizontalAdd rateMapExpanded];
            shiftedColCenter = colCenter + horizontalShift; 
        end
            
        if (yCenter > length(rateMapExpanded) / 2)
            shiftedRateMap = [shiftedRateMap; verticalAdd];
            shiftedRowCenter = rowCenter;
        else
            shiftedRateMap = [verticalAdd; shiftedRateMap];
            shiftedRowCenter = rowCenter + verticalShift; 
        end
        
        % Excludes the central peak from the grid score calculation 
        fieldWithoutCenter = shiftedRateMap;
        fieldWithoutCenter(shiftedRowCenter, shiftedColCenter) = 0;
        gridScore = calculate_grid_score(fieldWithoutCenter);
        if (gridScore > maxGridScore)
            maxGridScore = gridScore;
            maxShiftedField = shiftedRateMap;
        end
    end 
end

end
