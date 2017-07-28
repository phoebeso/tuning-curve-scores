clear all; clc

% load borderdata.mat
% rateMap = borderdata;
load simdata.mat

for i = 1:4
    rateMap = zeros(25);
    rateMap(1:15,1) = 1; 
    rateMap(1:15,2) = 0.5;
    rateMap(2:13,3) = 0.25;
    rateMap(7:14, 24:25) = 1;

%     rateMap = simdata{i};

%     rateMap = ones(25);
    
%     rateMap = zeros(25);
% %     rateMap(:,1) = 1;
%     rateMap(13,13) = 1;

    borderScore = calculate_border_score(rateMap);

    figure(i)
    imagesc(rateMap); colorbar
    xlabel(['Border Score: ' num2str(borderScore)])
end
