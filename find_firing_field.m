function [x_center, y_center, radius,rate_map] = find_firing_field(rate_map)
    peakrate = max(rate_map(:));
    threshold = 0.8 * peakrate;
    rate_map(rate_map <= threshold | isnan(rate_map)) = 0;
    rate_map(rate_map > threshold) = 1;
    rate_map = bwlabel(rate_map);
    
    % find groups of 6 bins???
    coordinates = [];
    for i = 1:max(rate_map(:))
        [row, col] = find(rate_map == i);
        if length(row) < 6
            rate_map(rate_map == i) = 0;
        else
            rate_map(rate_map == i) = 1;
            coordinates = vertcat(coordinates,[row col]);
        end
    end
    
    for i = 1:length(rate_map)
        for j = 1:length(rate_map)
            
        end
    end
    
    x_center = 0;
    y_center = 0;
    radius = 0;

end

