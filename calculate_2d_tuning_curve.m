function [tuningCurve, smoothTuningCurve] = calculate_2d_tuning_curve(variable_x,variable_y,fr,numBin,minVal,maxVal)
% Calculates a 2d spatial rate map 
% Assumes that the 2d environment is a square box, and that the
% variable is recorded along the x- and y-axes

%% define the axes and initialize variables

xAxis = linspace(minVal,maxVal,numBin+1);
yAxis = linspace(minVal,maxVal,numBin+1);

% initialize 
tuningCurve = zeros(numBin,numBin);

%% fill out the tuning curve

% find the mean firing rate in each position bin
for i  = 1:numBin
    start_x = xAxis(i); stop_x = xAxis(i+1);
    % find the times the animal was in the bin
    if i == numBin
        x_ind = find(variable_x >= start_x & variable_x <= stop_x);
    else
        x_ind = find(variable_x >= start_x & variable_x < stop_x);
    end
    
    for j = 1:numBin
        
        start_y = yAxis(j); stop_y = yAxis(j+1);
        
        if j == numBin
            y_ind = find(variable_y >= start_y & variable_y <= stop_y);
        else
            y_ind = find(variable_y >= start_y & variable_y < stop_y);
        end
        
        ind = intersect(x_ind,y_ind);
        
        % fill in rate map
        tuningCurve(numBin+1 - j,i) = mean(fr(ind));
    end
end

%% smooth the tuning curve

% fill in the NaNs with neigboring values
nan_ind = find(isnan(tuningCurve));
[j,i] = ind2sub(size(tuningCurve),nan_ind);
nan_num= numel(nan_ind);

% fill in the NaNs with neigboring values
for n = 1:nan_num
    ind_i = i(n); ind_j = j(n);
    
    right = tuningCurve(ind_j,min(ind_i+1,numBin));
    left = tuningCurve(ind_j,max(ind_i-1,1));
    down = tuningCurve(min(ind_j+1,numBin),ind_i);
    up = tuningCurve(max(ind_j-1,1),ind_i);
    
    ru = tuningCurve(max(ind_j-1,1),min(ind_i+1,numBin));
    lu = tuningCurve(max(ind_j-1,1),max(ind_i-1,1));
    ld = tuningCurve(min(ind_j+1,numBin),max(ind_i-1,1));
    rd = tuningCurve(max(ind_j-1,1),min(ind_i+1,numBin));
    
    tuningCurve(ind_j,ind_i) = nanmean([left right up down lu ru rd ld]);
    
end

% smooth with Gaussian
H = fspecial('gaussian'); % using default values - size=[3 3] and sigma=0.5
smoothTuningCurve = imfilter(tuningCurve,H);

return

