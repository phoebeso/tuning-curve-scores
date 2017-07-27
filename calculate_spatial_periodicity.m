function [rotations, correlations, circularMatrix] = calculate_spatial_periodicity(autocorrelationMatrix)
% Given an autocorrelation matrix, crops a circular region and then
% calculates the spatial periodicity by rotation the autocorrelation matrix
% in steps of 6 degrees and computing the correlation

dim = size(autocorrelationMatrix);
yDim = dim(1); xDim = dim(2);

threshold = 0.1;

modifiedMatrix = autocorrelationMatrix;
modifiedMatrix(modifiedMatrix <= threshold | isnan(modifiedMatrix)) = 0;
modifiedMatrix(modifiedMatrix > threshold) = 1;
modifiedMatrix = bwlabel(modifiedMatrix);

centerValue = modifiedMatrix(ceil(yDim / 2), ceil(xDim / 2));
[rowCenter, colCenter] = find(modifiedMatrix == centerValue);
xCenter = ceil((max(colCenter) + min(colCenter)) / 2);
yCenter = ceil((max(rowCenter) + min(rowCenter)) / 2);
center = [yCenter xCenter];

[row, col] = find(modifiedMatrix ~= 0 & modifiedMatrix ~= centerValue);
coordinates = [row col];
coordinates(:,3) = arrayfun(@(x,y) modifiedMatrix(x,y), row, col); % coordinate values in modifiedMatrix

% Compute radius of circular mask
% Radius defined as shortest distance of center to another field + the
% longest distance across any field

% Computes shortest Euclidean distance from center point 
distances = sqrt(sum(bsxfun(@minus, coordinates, center).^2,2));
minDistance = min(distances);

% Computes longest distance between any two points in any field 
nFields = max(modifiedMatrix(:));
maxLongestDistance = 0; 
for i = 1:nFields
    [row, col] = find(modifiedMatrix == i);
    longestDistance = calculate_longest_distance(col, row);
    maxLongestDistance = max(maxLongestDistance, longestDistance);
end

radius = minDistance + maxLongestDistance;

[xMask,yMask] = meshgrid(-(xCenter-1):(xDim-xCenter),-(yCenter-1):(yDim-yCenter));
mask = ((xMask.^2 + yMask.^2) <= radius^2);
circularMatrix = autocorrelationMatrix .* mask;
 
% Concatanates zero vectors horizontally and vertically to center the
% rate map around the field of interest for later crosscorrelation 
% calculations
% Additionally shifts the central field's rows and cols for
% identification purposes for calculations
horizontalShift = abs(xCenter - xDim + xCenter);
verticalShift = abs(yCenter - yDim + yCenter);
horizontalAdd = zeros(yDim, horizontalShift);
verticalAdd = zeros(verticalShift, xDim + horizontalShift);
if (xCenter > xDim / 2)
    shiftedCircularMatrix = [circularMatrix horizontalAdd];
    shiftedColCenter = colCenter;
else
    shiftedCircularMatrix = [horizontalAdd circularMatrix];
    shiftedColCenter = colCenter + horizontalShift; 
end

if (yCenter > yDim / 2)
    shiftedCircularMatrix = [shiftedCircularMatrix; verticalAdd];
    shiftedRowCenter = rowCenter;
else
    shiftedCircularMatrix = [verticalAdd; shiftedCircularMatrix];
    shiftedRowCenter = rowCenter + verticalShift; 
end

rotations = 0:6:360;
correlations = nan(61,1);
for i = 1:61
    correlations(i) = calculate_correlation(shiftedCircularMatrix, 6*(i-1));
end

end
