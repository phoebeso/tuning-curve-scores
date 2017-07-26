% Border scores are calculated from two variables, CM and DM. First, firing
% fields were defined as group of adjacent pixels with a firing rate larger
% than 20% of the peak firing rate of the rate map. For each field, the
% proportion of pixels directly adjacent to a wall and also part of the
% field was calculated. CM was defined as the maximum proportion obtained
% over all fields. DM was the average shortest distance to a wall for all
% pixels part of a firing field, weighted by the firing rate in each pixel.
% DM was then normalized by the longest shortest distance to a wall of any
% firing pixel in the rate map. The border score was then defined as 
% (CM - DM)/(CM + DM)
% Border scores range from +1 for a rate map with one infinitely thin
% firing field along an entire wall, to -1 for a rate map with an
% infinitely small firing field in the center of the map. 

clear all; clc

% load borderdata.mat
% rateMap = borderdata;
load simdata.mat

rateMap = simdata{1};
rateMap(:,7:18) = 0;

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
maxProportion = 0;
for i = 1:nFields
    % Determines number of adjacent pixels to a wall in a firing field
    [fieldRow, fieldCol] = find(modifiedRateMap == i);
    nFieldPixels = length(fieldRow);
    % Fields defined as needing to have more than 6 pixels 
    if (nFieldPixels <= 6)
        modifiedRateMap(fieldRow, fieldCol) = 0;
        continue;
    end
    
    fieldWidth = max(fieldCol) - min(fieldCol) + 1;
    fieldHeight = max(fieldRow) - min(fieldRow) + 1;
    % pixels along vertical wall
    if (fieldHeight > fieldWidth)
        if (mean(fieldCol) < length(rateMap) / 2)
            nAdjacentPixels = length(find(fieldCol == 1)); % pixels along left wall
        else
            nAdjacentPixels = length(find(fieldCol == length(rateMap))); % pixels along right wall
        end
    % pixels along horizontal wall 
    else
        if (mean(fieldRow) < length(rateMap) / 2)
            nAdjacentPixels = length(find(fieldRow == 1)); % pixels along top wall
        else
            nAdjacentPixels = length(find(fieldRow == length(rateMap))); % pixels along bottom wall
        end
    end
    
    proportion = nAdjacentPixels / nFieldPixels;
    maxProportion = max(maxProportion, proportion); 
end

CM = maxProportion; 

% Calculates DM
[row, col] = find(modifiedRateMap ~= 0);
nPixels = length(row);
weightedTotalDistance = 0;
maxDistance = 0;
for i = 1:nPixels
    % Calculates distance of pixel from all four walls
    distances = [abs(row(i) - 1) abs(row(i) - length(rateMap)) abs(col(i) - 1) abs(col(i) - length(rateMap))];
    minDistance = min(distances);
    rate = rateMap(row(i), col(i));
    weight = rate / maxRate; 
    weightedDistance = minDistance * weight;
    weightedTotalDistance = weightedTotalDistance + weightedDistance;
    % Largest smallest distance of any pixel from a wall 
    maxDistance = max(maxDistance, minDistance);
end

weightedAverageDistance = weightedTotalDistance / nPixels;

DM = weightedAverageDistance / maxDistance;
if (isnan(DM))
    DM = 0;
end

borderScore = (CM - DM)/(CM + DM); 
