function [collapsePartitionData] = analyze_periodicity(rotations,correlations)
% Correlations = correlation of the autocorrelation matrix every 6 degrees

rotations = rotations(1:60);
correlations = correlations(1:60);

cellIdx = 1;
collapsePartitionData = cell(8,2);
for nPartitions = 3:10
    periods = linspace(0,360,nPartitions+1);
    [~, ~, bin] = histcounts(rotations, periods);
    partitions = zeros(nPartitions,floor(360/6/nPartitions));
    for i = 1:nPartitions
        partitionData = correlations(bin == i)';
        partitions(i, :) = partitionData(1:floor(360/6/nPartitions));
    end
    
    sumPartitions = sum(partitions);
    
    collapsePartitionData{cellIdx, 1} = nPartitions; % number of partitions 
    collapsePartitionData{cellIdx, 2} = sumPartitions; % collapsed partitons data 
    
    cellIdx = cellIdx + 1;
end

end

