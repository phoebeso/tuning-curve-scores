function [crosscorrelation] = calculate_crosscorrelation(graph, angle)
% Rotates rate map by specified angle 
rotated_graph = imrotate(graph,angle,'crop');
rotated_graph(rotated_graph == 0) = NaN;
    
% Pixels in the rate map for which rate was estimated for both the original
% and rotated rate map
n = sum(sum(~isnan(graph)));
    
% Calculates the spatial crosscorrelation of the original and rotated rate
% map
[sum1, sum2, sum3, sum4, sum5] = deal(0);
for i = 1:length(graph)
    for j = 1:length(graph)
        if ~isnan(graph(i,j)) && ~isnan(rotated_graph(i,j))
            sum1 = sum1 + (graph(i,j) * rotated_graph(i,j));
            sum2 = sum2 + graph(i,j);
            sum3 = sum3 + rotated_graph(i,j);
            sum4 = sum4 + graph(i,j)^2;
            sum5 = sum5 + rotated_graph(i,j)^2;
        else
            n = n - 1;
        end
    end
end
    
numerator = (n * sum1) - (sum2 * sum3);
denominator = sqrt((n * sum4) - sum2^2) * sqrt((n * sum5) - sum3^2);
crosscorrelation = numerator / denominator;

end
