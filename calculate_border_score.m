function [borderScore] = calculate_border_score(rateMap)
% We identified border cells by computing, for each cell, the difference between 
% the maximal length of a wall touching on a single firing field and the average 
% distance of the firing locations from the nearest wall, divided by the sum of those values
% http://www.nature.com/neuro/journal/v13/n8/full/nn.2602.html

% OR

% We identified such cells by computing, for each cell, the difference between 
% the maximal length of a single wall touching on a single firing field and 
% the average distance of this field from the wall, divided by the sum of those values
% http://www.sciencedirect.com/science/article/pii/S0896627314001123

maxRate = max(rateMap(:));
threshold = 0.2 * maxRate;

% Converts all position bins with firing rates below the threshold to zero
% and labels connected components with firing rates above the threshold
modifiedRateMap = rateMap;
modifiedRateMap(modifiedRateMap <= threshold | isnan(modifiedRateMap)) = 0;
modifiedRateMap(modifiedRateMap > threshold) = 1;
modifiedRateMap = bwlabel(modifiedRateMap);

% Calculates CM
nFields = max(modifiedRateMap(:));
maxDistance = 0;
longestField = 0;
maxIdx = 0;
for i = 1:nFields
    % Determines number of adjacent pixels to a wall in a firing field
    [fieldRow, fieldCol] = find(modifiedRateMap == i);
    
    nFieldPixels = length(fieldRow);
    % Fields defined as needing to have more than 6 pixels 
    if (nFieldPixels <= 6)
%         modifiedRateMap(modifiedRateMap == i) = 0;
        continue;
    end
    
    % Calculates length of pixels along each wall 
    [rowLeftPixels, ~] = find(fieldCol == 1); % pixels along vertical left wall
    leftDistance = abs(max(rowLeftPixels) - min(rowLeftPixels));
    
    [rowRightPixels, ~] = find(fieldCol == length(rateMap)); % pixels along vertical right wall
    rightDistance = abs(max(rowRightPixels) - min(rowRightPixels));
    
    [~, colTopPixels] = find(fieldRow == 1); % pixels along horizontal top wall
    topDistance = abs(max(colTopPixels) - min(colTopPixels));
    
    [~, colBottomPixels] = find(fieldRow == length(rateMap)); % pixels along horizontal bottom wall
    bottomDistance = abs(max(colBottomPixels) - min(colBottomPixels));

    [fieldMaxDistance, idx] = max([leftDistance rightDistance topDistance bottomDistance]);
    
    if fieldMaxDistance > maxDistance
        longestField = i;
        maxIdx = idx; 
        maxDistance = fieldMaxDistance;
    end
end

CM = maxDistance;

% Calculates DM
[row, col] = find(modifiedRateMap == longestField);
% [row, col] = find(modifiedRateMap ~= 0); 
nPixels = length(row);
totalDistance = 0;
for i = 1:nPixels
    if maxIdx == 1 || maxIdx == 3
        wall = 1;
    else
        wall = length(rateMap);
    end
    
    if maxIdx == 1 || maxIdx == 2
        distance = abs(col(i) - wall);
    else
        distance = abs(row(i) - wall);
    end
    
    totalDistance = totalDistance + distance;
end

DM = totalDistance / nPixels;
if (isnan(DM))
    DM = 0;
end

borderScore = (CM - DM)/(CM + DM);

end
