clear all; clc

% load borderdata.mat
% rateMap = borderdata;
load simdata.mat

for i = 1:4
    rateMap = zeros(25);
    rateMap(1:15,1) = 1; 
    rateMap(1:15,2) = 0.5;
    rateMap(2:13,3) = 0.25;
% %     rateMap(1:15,1:3) = ones(15,3);
% %     rateMap(18:19, 16:17) = ones(2,2);
    borderScore = calculate_border_score(rateMap);
%     rateMap = simdata{i};
%     borderScore = calculate_border_score(rateMap);

    figure(i)
    imagesc(rateMap); colorbar
    xlabel(['Border Score: ' num2str(borderScore)])
end
