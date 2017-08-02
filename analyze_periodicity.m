function [collapsePartitionData] = analyze_periodicity(rotations,correlations)
% Partitions the sinusoidal correlation data into n even partitions. Collapses the
% partitions by calculating the summation of all n partitions. Data is used
% for periodicity analysis.

rotations = rotations(1:60);
correlations = correlations(1:60);

% Partitions the sinusoidal correlation data in 3-10 even partitions and
% collapses the data
cellIdx = 1;
collapsePartitionData = cell(8,3);
for nPartitions = 3:10
    periods = linspace(0,360,nPartitions+1);
    [~, ~, bin] = histcounts(rotations, periods);
    
    % Gets all n partitions
    partitions = zeros(nPartitions,floor(360/6/nPartitions));
    for i = 1:nPartitions
        partitionData = correlations(bin == i)';
        partitions(i, :) = partitionData(1:floor(360/6/nPartitions));
    end
    
    % Collapses partitions
    sumPartitions = sum(partitions);
    
    collapsePartitionData{cellIdx, 1} = nPartitions; % number of partitions
    collapsePartitionData{cellIdx, 2} = max(sumPartitions); % max value of collapsed data 
    collapsePartitionData{cellIdx, 3} = sumPartitions; % collapsed partitons data 
    
    cellIdx = cellIdx + 1;
end

end

