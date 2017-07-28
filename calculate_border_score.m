function [borderScore] = calculate_border_score(rateMap)
% Border scores are calculated from two variables, CM and DM. CM was
% defined as the maximal length of a single wall touching on a single
% firing field. DM was defined as the average minimum distance of all 
% pixels in a field to any wall. The border score was then defined as 
% (CM - DM)/(CM + DM)

maxRate = max(rateMap(:));
threshold = 0.2 * maxRate;

% Converts all position bins with firing rates below the threshold to zero
% and labels connected components with firing rates above the threshold
modifiedRateMap = rateMap;
modifiedRateMap(modifiedRateMap <= threshold | isnan(modifiedRateMap)) = 0;
modifiedRateMap(modifiedRateMap > threshold) = 1;
modifiedRateMap = bwlabel(modifiedRateMap);
nFields = max(modifiedRateMap(:));

% Calculates CM. CM is the largest proportion of pixels directly next to a
% wall in a single firing field
% maxProportion = 0;
maxDistance = 0;
for i = 1:nFields
    % Determines number of adjacent pixels to a wall in a firing field
    [fieldRow, fieldCol] = find(modifiedRateMap == i);
    
    nFieldPixels = length(fieldRow);
    % Fields defined as needing to have more than 6 pixels 
    if (nFieldPixels <= 6)
        modifiedRateMap(modifiedRateMap == i) = 0;
        continue;
    end
    
    % Calculates number of pixels along each wall and identifies the wall
    % with the highest number of adjacent pixels  
    [rowLeftPixels, ~] = find(fieldCol == 1); % pixels along vertical left wall
    leftDistance = abs(max(rowLeftPixels) - min(rowLeftPixels));
    
    [rowRightPixels, ~] = find(fieldCol == length(rateMap)); % pixels along vertical right wall
    rightDistance = abs(max(rowRightPixels) - min(rowRightPixels));
    
    [~, colTopPixels] = find(fieldRow == 1); % pixels along horizontal top wall
    topDistance = abs(max(colTopPixels) - min(colTopPixels));
    
    [~, colBottomPixels] = find(fieldRow == length(rateMap)); % pixels along horizontal bottom wall
    bottomDistance = abs(max(colBottomPixels) - min(colBottomPixels));

    maxDistance = max([maxDistance leftDistance rightDistance topDistance bottomDistance]);
end

CM = maxDistance;

% Calculates DM
[row, col] = find(modifiedRateMap ~= 0);
nPixels = length(row);
totalDistance = 0;
for i = 1:nPixels
    distances = [abs(row(i) - 1) abs(row(i) - length(rateMap)) abs(col(i) - 1) abs(col(i) - length(rateMap))];
    minDistance = min(distances);
    
    totalDistance = totalDistance + minDistance;
end

DM = totalDistance / nPixels;
if (isnan(DM))
    DM = 0;
end

borderScore = (CM - DM)/(CM + DM); 

end

