function [rotations, correlations, matrixWithoutCenter] = calculate_spatial_periodicity(autocorrelationMatrix)
% Given an autocorrelation matrix, crops a circular region and then
% calculates the spatial periodicity by rotation the autocorrelation matrix
% in steps of 6 degrees and computing the correlation

dim = size(autocorrelationMatrix);
yDim = dim(1); xDim = dim(2);

threshold = 0.1;

modifiedMatrix = autocorrelationMatrix;
modifiedMatrix(modifiedMatrix <= threshold | isnan(modifiedMatrix)) = 0;
modifiedMatrix(modifiedMatrix > threshold) = 1;
modifiedMatrix = bwlabel(modifiedMatrix); % IDs the fields with a value 

% Determines and identifies the center of the autocorrelation matrix
yMatrixCenter = ceil(yDim / 2); xMatrixCenter = ceil(xDim / 2);
if (modifiedMatrix(yMatrixCenter, xMatrixCenter) ~= 0)
    centerValue = modifiedMatrix(ceil(yDim / 2), ceil(xDim / 2));
else
    % Coordinates of all fields in the autocorrelation matrix
    [rowAllFields, colAllFields] = find(modifiedMatrix ~= 0);
    % Distances of all coordinates from the center of the matrix 
    distanceFromCenter = sqrt(sum(bsxfun(@minus, [rowAllFields colAllFields], [yMatrixCenter xMatrixCenter]).^2,2));
    [~, idx1] = min(distanceFromCenter);
    centerValue = modifiedMatrix(rowAllFields(idx1(1)), colAllFields(idx1(1)));
end

% Gets center of central field 
[rowCenterField, colCenterField] = find(modifiedMatrix == centerValue);
xCenterField = ceil((max(colCenterField) + min(colCenterField)) / 2);
yCenterField = ceil((max(rowCenterField) + min(rowCenterField)) / 2);
center = [yCenterField xCenterField];

% Creates a matrix of all coordinates in fields besides the central field.
% First col is the row of the coordinate, second col is the col, third
% col is the ID value of the coordinate in the modifiedMatrix (IDs which 
% field the coordinate belongs to), and fourth col is the distance of the 
% coordinate to the center of the central field
[rowFields, colFields] = find(modifiedMatrix ~= 0 & modifiedMatrix ~= centerValue);
coordinates = [rowFields colFields];
coordinates(:,3) = arrayfun(@(x,y) modifiedMatrix(x,y), rowFields, colFields);
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
    fieldValues = sortedCoordinates(sortedIdx2(1:6),3);
else
    fieldValues = sortedCoordinates(sortedIdx2, 3);
end

[rowCircularFields, colCircularFields] = find(ismember(modifiedMatrix,fieldValues));
circularFieldsCoord = [rowCircularFields colCircularFields];
distances = sqrt(sum(bsxfun(@minus, circularFieldsCoord, center).^2,2));

radius = max(distances);

[xMask,yMask] = meshgrid(-(xCenterField-1):(xDim-xCenterField),-(yCenterField-1):(yDim-yCenterField));
mask = ((xMask.^2 + yMask.^2) <= radius^2);
mask = double(mask);
mask(mask == 0) = NaN;
circularMatrix = autocorrelationMatrix .* mask;
 
% Concatanates zero vectors horizontally and vertically to center the
% rate map around the field of interest for later crosscorrelation 
% calculations
% Additionally shifts the central field's rows and cols for
% identification purposes for calculations
horizontalShift = abs(xCenterField - xDim + xCenterField);
verticalShift = abs(yCenterField - yDim + yCenterField);
horizontalAdd = nan(yDim, horizontalShift);
verticalAdd = nan(verticalShift, xDim + horizontalShift);
if (xCenterField > xDim / 2)
    shiftedCircularMatrix = [circularMatrix horizontalAdd];
    shiftedColCenter = colCenterField;
else
    shiftedCircularMatrix = [horizontalAdd circularMatrix];
    shiftedColCenter = colCenterField + horizontalShift; 
end

if (yCenterField > yDim / 2)
    shiftedCircularMatrix = [shiftedCircularMatrix; verticalAdd];
    shiftedRowCenter = rowCenterField;
else
    shiftedCircularMatrix = [verticalAdd; shiftedCircularMatrix];
    shiftedRowCenter = rowCenterField + verticalShift; 
end

matrixWithoutCenter = shiftedCircularMatrix;
matrixWithoutCenter(shiftedRowCenter, shiftedColCenter) = NaN;
rotations = 0:6:360;
correlations = nan(61,1);
for i = 1:61
    correlations(i) = calculate_correlation(matrixWithoutCenter, 6*(i-1));
end

end
