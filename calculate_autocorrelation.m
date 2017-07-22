function [autocorrelation] = calculate_autocorrelation(rateMap, angle)
% Calculates the correlation between a graph and the graph rotated by a
% specified angle

% Rotates rate map by specified angle 
rotatedGraph = imrotate(rateMap,angle,'crop');
rotatedGraph(rotatedGraph == 0) = NaN;
    
% Number of pixels in the rate map for which rate was estimated for both 
% the original and rotated rate map
n = sum(~isnan(rateMap(:)));
    
% Calculates the spatial correlation of the original and rotated rate map
[sum1,sum2,sum3,sum4,sum5] = deal(0);
dimensions = size(rateMap);
for i = 1:dimensions(1)
    for j = 1:dimensions(2)
        if (~isnan(rateMap(i,j))) && (~isnan(rotatedGraph(i,j)))
            sum1 = sum1 + (rateMap(i,j) * rotatedGraph(i,j));
            sum2 = sum2 + rateMap(i,j);
            sum3 = sum3 + rotatedGraph(i,j);
            sum4 = sum4 + rateMap(i,j)^2;
            sum5 = sum5 + rotatedGraph(i,j)^2;
        else
            n = n - 1;
        end
    end
end
    
numerator = (n * sum1) - (sum2 * sum3);
denominator = sqrt((n * sum4) - sum2^2) * sqrt((n * sum5) - sum3^2);
autocorrelation = numerator / denominator;

end
