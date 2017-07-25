function [autocorrelogram] = calculate_autocorrelogram(graph)
% Generates the autocorrelogram by averaging the rates of the graph rotated
% every 60 degrees
graph(isnan(graph)) = 0;
autocorrelogram = graph;
nRotations = 5;
for i = 1:nRotations
    angle = 60 * i; 
    rotatedGraph = imrotate(graph,angle,'crop');
    rotatedGraph(isnan(rotatedGraph)) = 0;
    autocorrelogram = autocorrelogram + rotatedGraph; 
end

autocorrelogram = autocorrelogram ./ (nRotations + 1); 

end
