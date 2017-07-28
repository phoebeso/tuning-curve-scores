function [correlation] = calculate_correlation(matrix, angle)
% Calculates the correlation between a graph and the graph rotated by a
% specified angle

% How much the rate map was expanded by (look at find_central_field line 17 code)
% expand = 3;

% Rotates rate map. Also resolves issue of a position having a value of 0
minValue = min(matrix(:));
tempMatrix = matrix - minValue + 1; % Sets minimum to 1 
rotatedMatrix = imrotate(tempMatrix, angle, 'crop');
rotatedMatrix(rotatedMatrix == 0) = NaN;
rotatedMatrix = rotatedMatrix + minValue - 1; % Readjusts values 

% Number of pixels in the rate map for which rate was estimated for both 
% the original and rotated rate map
n = 0;

% Calculates the correlation of the original and rotated rate map
dimensions = size(matrix);
[sum1, sum2, sum3, sum4, sum5] = deal(0);
for i = 1:dimensions(1)
    for j = 1:dimensions(2)
        if (~isnan(matrix(i,j)) && ~isnan(rotatedMatrix(i,j)))
            sum1 = sum1 + (matrix(i,j) * rotatedMatrix(i,j));
            sum2 = sum2 + matrix(i,j);
            sum3 = sum3 + rotatedMatrix(i,j);
            sum4 = sum4 + matrix(i,j)^2;
            sum5 = sum5 + rotatedMatrix(i,j)^2;
            n = n + 1;
        end
    end
end

% Only estimate autocorrelation if n > 20 * expand pixels
% if (n < 20 * expand)
if (n < 20)
    correlation = NaN;
else 
    numerator = (n * sum1) - (sum2 * sum3);
    denominator = sqrt((n * sum4) - sum2^2) * sqrt((n * sum5) - sum3^2);
    correlation = numerator / denominator;
end

end
