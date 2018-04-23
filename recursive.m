function [pins, addedPins] = recursive(pins, index, max_rows, max_col, addedPins)
%Recursive function that solves for the cluster
    
    %Given the index of the broken pin
    %Find the broken neighbors
    i = ceil(index/max_rows);
    j = mod(index-1,max_col)+1;
    
    for m = 1:3
        for n = 1:3
            %Find the neighbor's value
            iNB = i + (m-2);
            jNB = j + (n-2);
            
            index = (iNB-1)*max_rows+jNB;
            Loc = find(index == pins(1,:));
            
            %If the the new node isn't already pinned
            if ~isempty(Loc) && pins(8,Loc) ~= 2
                pins(8,Loc) = 2; %Save as temporary pin
                addedPins = [addedPins; index];
                [pins addedPins] = recursive(pins, index, max_rows, max_col, addedPins);
            end
        end
    end
end

