function [theta, r] = adjusted_polar_graph(nPartitions, partitionData)
% Plots a polar bar graph of the collapsed partition data and adjusts the
% graph's angle labels to match the correct period

% Expands partitionData so data will be plotted as a bar instead of a point
r = [repmat(partitionData, 6*nPartitions, 1); zeros(1,length(partitionData))];
r = [0 reshape(r, 1, [])];
        
theta = linspace(0, 360, length(r));
theta = deg2rad(theta);

polar(theta,r);

% Deletes default angle labels on polar graph and replaces them with
% adjusted values 
hiddenText = findall(gca,'type','text');
angles = 0:30:330;
hObjToDelete = zeros(length(angles)-4,1);
angleIncrement = 360/nPartitions/4;
i = 0;
for ang = angles
    hObj = findall(hiddenText,'string',num2str(ang));
    switch ang
        case 0
            set(hObj,'string',sprintf('%.1f%c', 0, char(176)));
        case 90
            set(hObj,'string',sprintf('%.1f%c', angleIncrement, char(176)));
        case 180
            set(hObj,'string',sprintf('%.1f%c', angleIncrement*2, char(176)));
        case 270
            set(hObj,'string',sprintf('%.1f%c', angleIncrement*3, char(176)));
        otherwise
            i = i + 1;
            hObjToDelete(i) = hObj;
    end
end
delete(hObjToDelete);

title(sprintf('%.1f%c Period', 360/nPartitions, char(176)))

end

