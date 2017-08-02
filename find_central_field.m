function [maxShiftedRateMap, maxGridScore] = find_central_field(rateMap)
% Given a rate map, determines the center field to rotate the map around
% from multiple firing fields that maximizes the grid score 

% Threshold firing rate for a position bin to be considered part of a field
maxRate = max(rateMap(:));
threshold = 0.2 * maxRate;
 
% Converts all position bins with firing rates below the threshold to zero
% and labels connected components with firing rates above the threshold
modifiedRateMap = rateMap;
modifiedRateMap(modifiedRateMap <= threshold | isnan(modifiedRateMap)) = 0;
modifiedRateMap(modifiedRateMap > threshold) = 1;
modifiedRateMap = bwlabel(modifiedRateMap);
        
% Increases number of pixels per time bin to increase image/circle
% resolution
expand = 3;
modifiedRateMapExpanded = kron(modifiedRateMap, ones(expand)); % Map of peaks 
rateMapExpanded = kron(rateMap, ones(expand)); % Original map of all firing rates
 
% Loops through all fields and centers the rate map around the field that
% yields the highest grid score
maxGridScore = -inf;
maxShiftedRateMap = rateMapExpanded;
for i = 1:max(modifiedRateMap(:))
    % Only considers connected components with greater than 6*expand
    % position bins as fields
    % 6 arbitrarily chosen, but normally a firing field was defined as
    % being a contiguous region of at least 225 cm^2
    fieldSize = length(find(modifiedRateMapExpanded == i));
    if (fieldSize > 6 * expand)
        % Gets coordinates of the field to calculate the center of
        % rotation
        [rowCenter, colCenter] = find(modifiedRateMapExpanded == i);
        xCenter = ceil((max(colCenter) + min(colCenter)) / 2);
        yCenter = ceil((max(rowCenter) + min(rowCenter)) / 2);
           
        % Concatanates zero vectors horizontally and vertically to center the
        % rate map around the field of interest for later crosscorrelation 
        % calculations
        % Additionally shifts the central field's rows and cols for
        % identification purposes for calculations
        dim = size(rateMapExpanded);
        yDim = dim(1); xDim = dim(2);
        horizontalShift = abs(xCenter - xDim + xCenter);
        verticalShift = abs(yCenter - yDim + yCenter);
        horizontalAdd = nan(yDim, horizontalShift);
        verticalAdd = nan(verticalShift, xDim + horizontalShift);
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
        
        % Excludes the central field from the grid score calculation 
        fieldWithoutCenter = shiftedRateMap;
        centerIdx = sub2ind(size(fieldWithoutCenter), shiftedRowCenter, shiftedColCenter);
        fieldWithoutCenter(centerIdx) = NaN;
        gridScore = calculate_grid_score(fieldWithoutCenter);
        if (gridScore > maxGridScore)
            maxGridScore = gridScore;
            maxShiftedRateMap = shiftedRateMap;
        end
    end
end

end
