function [rotations, correlations, matrixWithoutCenter] = calculate_spatial_periodicity(autocorrelationMatrix)
% Given an autocorrelation matrix, crops a circular region and then
% calculates the spatial periodicity by rotating the autocorrelation matrix
% in steps of 6 degrees and computing the correlation

% Locates peaks in the autocorrelation matrix and groups them together  
threshold = 0.1;
modifiedMatrix = autocorrelationMatrix;
modifiedMatrix(modifiedMatrix <= threshold | isnan(modifiedMatrix)) = 0;
modifiedMatrix(modifiedMatrix > threshold) = 1;
modifiedMatrix = bwlabel(modifiedMatrix); % IDs the pekas with a value

% Enhances matrices for smoother circle
expand = 3;

modifiedMatrix = kron(modifiedMatrix, ones(expand));
autocorrelationMatrix = kron(autocorrelationMatrix, ones(expand)); 

dim = size(autocorrelationMatrix);
yDim = dim(1); xDim = dim(2);

% Determines and identifies the center of the autocorrelation matrix
yMatrixCenter = ceil(yDim / 2); xMatrixCenter = ceil(xDim / 2);
if (modifiedMatrix(yMatrixCenter, xMatrixCenter) ~= 0)
    centerId = modifiedMatrix(ceil(yDim / 2), ceil(xDim / 2));
else
    % Coordinates of all fields in the autocorrelation matrix
    [rowAllPeaks, colAllPeaks] = find(modifiedMatrix ~= 0);
    % Distances of all coordinates from the center of the matrix 
    distanceFromCenter = sqrt(sum(bsxfun(@minus, [rowAllPeaks colAllPeaks], [yMatrixCenter xMatrixCenter]).^2,2));
    [~, idx1] = min(distanceFromCenter);
    centerId = modifiedMatrix(rowAllPeaks(idx1(1)), colAllPeaks(idx1(1)));
end

% Gets center of central peak 
[rowCenterPeak, colCenterPeak] = find(modifiedMatrix == centerId);
xCenterPeak = ceil((max(colCenterPeak) + min(colCenterPeak)) / 2);
yCenterPeak = ceil((max(rowCenterPeak) + min(rowCenterPeak)) / 2);
center = [yCenterPeak xCenterPeak];

% Creates a matrix of all coordinates in fields besides the central field.
% First col is the row of the coordinate, second col is the col, third
% col is the ID value of the coordinate in the modifiedMatrix (IDs which 
% field the coordinate belongs to), and fourth col is the distance of the 
% coordinate to the center of the central field
[rowPeaks, colPeaks] = find(modifiedMatrix ~= 0 & modifiedMatrix ~= centerId);
coordinates = [rowPeaks colPeaks];
coordinates(:,3) = arrayfun(@(x,y) modifiedMatrix(x,y), rowPeaks, colPeaks);
coordinates(:,4) = sqrt(sum(bsxfun(@minus, coordinates(:,1:2), center).^2,2));

% Sorts coordinates by distance to center in order to identify the 6 cloest
% fields to the center field 
[~, idx2] = sort(coordinates(:,4));
sortedCoordinates = coordinates(idx2,:);
% Determines 6 closest fields by looking at unique field ID values 
[~, idx3, ~] = unique(sortedCoordinates(:,3), 'first');
sortedIdx2 = sort(idx3);
% Gets field ID values of 6 closest values
if (length(sortedIdx2) >= 6)
    PeakIds = sortedCoordinates(sortedIdx2(1:6),3);
else
    PeakIds = sortedCoordinates(sortedIdx2, 3);
end

% Gets coordinats of 6 closest fields and calculates the radius of the
% circular area as the farthest distance from the center to any coordinate
% in the 6 fields
[rowCircularPeaks, colCircularPeaks] = find(ismember(modifiedMatrix,PeakIds));
circularPeaksCoord = [rowCircularPeaks colCircularPeaks];
distances = sqrt(sum(bsxfun(@minus, circularPeaksCoord, center).^2,2));
radius = max(distances);

% Extracts circular area from autocorrelation matrix 
[xMask,yMask] = meshgrid(-(xCenterPeak-1):(xDim-xCenterPeak),-(yCenterPeak-1):(yDim-yCenterPeak));
mask = ((xMask.^2 + yMask.^2) <= radius^2);
% mask = double(mask);
% mask(mask == 0) = NaN;

% DISTANCE OF CENTER TO ANY OTHER POINT IN THE CENTRAL FIELD
centerDistances = sqrt(sum(bsxfun(@minus, [rowCenterPeak colCenterPeak], center).^2,2));
radius2 = max(centerDistances);
[xMask2,yMask2] = meshgrid(-(xCenterPeak-1):(xDim-xCenterPeak),-(yCenterPeak-1):(yDim-yCenterPeak));
mask2 = ((xMask2.^2 + yMask2.^2) > radius2^2);
% mask2 = double(mask2);

mask3 = mask & mask2;
mask3 = double(mask3);
mask3(mask3 == 0) = NaN;

circularMatrix = autocorrelationMatrix .* mask3;
 
% Concatanates nan vectors horizontally and vertically to center the
% circular area for later crosscorrelation calculations
% Additionally shifts the central field's rows and cols for
% identification purposes for calculations
horizontalShift = abs(xCenterPeak - xDim + xCenterPeak);
verticalShift = abs(yCenterPeak - yDim + yCenterPeak);
horizontalAdd = nan(yDim, horizontalShift);
verticalAdd = nan(verticalShift, xDim + horizontalShift);
if (xCenterPeak > xDim / 2)
    shiftedCircularMatrix = [circularMatrix horizontalAdd];
    shiftedColCenter = colCenterPeak;
else
    shiftedCircularMatrix = [horizontalAdd circularMatrix];
    shiftedColCenter = colCenterPeak + horizontalShift; 
end

if (yCenterPeak > yDim / 2)
    shiftedCircularMatrix = [shiftedCircularMatrix; verticalAdd];
    shiftedRowCenter = rowCenterPeak;
else
    shiftedCircularMatrix = [verticalAdd; shiftedCircularMatrix];
    shiftedRowCenter = rowCenterPeak + verticalShift; 
end

% Removes center field for correlation calculations
matrixWithoutCenter = shiftedCircularMatrix;
centerIdx = sub2ind(size(matrixWithoutCenter), shiftedRowCenter, shiftedColCenter);
matrixWithoutCenter(centerIdx) = NaN;

% Calculates the correlation every 6 degrees 
rotations = 0:6:360;
correlations = nan(61,1);
for i = 1:61
    correlations(i) = calculate_correlation(matrixWithoutCenter, 6*(i-1));
end

end
