function [CircCurr, MLVandBL, MLVnodes] = MLV(CircCurr, max_rows)
%SOURCE This function stores the index values of the MLV

    %Keep source value at 1 (set color in colormap)
    MLVandBL = CircCurr; %Store current expansion as the total output
    
    [i, j] = find(CircCurr); %Find the coordinates of the MLV nodes
    MLVSize = length(i); %Store the number of MLV nodes
    MLVnodes = zeros(1,MLVSize); 
    
    for h = 1:MLVSize
        index = indxx(i(h), j(h), max_rows);
        MLVnodes(1,h) = index; %Save the indexes of the MLV nodes
    end
end

