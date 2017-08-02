clear all; clc

load data_for_cell77.mat

posx = posx_c;
posy = posy_c;
nPosBins = 20;
boxSize = 100;

filter = gaussian(-4:4,2,0); filter = filter / sum(filter); 
dt = post(3) - post(2); firingRate = spiketrain / dt;
smoothFiringRate = conv(firingRate,filter,'same');

[rateMap] = compute_2d_tuning_curve(posx,posy,smoothFiringRate,nPosBins,0,boxSize);

load simdata.mat

for i = 1:1

    rateMap = simdata{1};

    [~, gridScore] = find_central_field(rateMap);

    matrix = calculate_correlation_matrix(rateMap);
    
    [rotations, correlations, circularMatrix] = calculate_spatial_periodicity(matrix);
    
    figure(i)
    subplot(2,2,1)
    imagesc(rateMap); colorbar
    title('Firing Rate Map')
    xlabel(['Grid Score: ' num2str(gridScore)])

    subplot(2,2,2)
    imagesc(matrix, [-1 1]); colorbar
    axis off
    title('Autocorrelation Matrix')

    subplot(2,2,3)
    imagesc(circularMatrix, [-1 1]); colorbar
    axis off
    title('Circular Autocorrelation Matrix')

    gridScore2 = calculate_grid_score(circularMatrix); 
    
	subplot(2,2,4)
    plot(rotations, correlations);
    xlim([0 360])
    ylim([-inf 1])
    xlabel({'Rotation (deg)'; sprintf('Grid Score: %f', gridScore2)})
    ylabel('Correlation')
    title('Periodicity')
    
    collapsePartitionData = analyze_periodicity(rotations, correlations);
    
    figure(i*2)
    for j = 1:length(collapsePartitionData)
        subplot(2,4,j)
        partitionData = collapsePartitionData{j,2};
        xr = [repmat(partitionData, 6*collapsePartitionData{j,1}, 1); zeros(1,length(partitionData))];
        xr = [0 reshape(xr, 1, [])];
        th = linspace(0, 359, length(xr));
        th = deg2rad(th);

        polar(th, xr);

        title(sprintf('%.1f%c Degree Period', 360/collapsePartitionData{j,1}, char(176)))
    end
    
end
