clear all; clc

load simdata.mat

for i = 1:4
%     rateMap = simdata{i};

    borderScore = calculate_border_score(rateMap);

    figure(i)
    imagesc(rateMap); colorbar
    xlabel(['Border Score: ' num2str(borderScore)])
end
