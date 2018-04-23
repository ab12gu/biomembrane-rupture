function [pins, LNBindex, MLVandBL] = ...
    cluster(pins, P, max_rows, max_col, MLVandBL, threshold, ycen, xcen,...
    maxRadius,  unpinned0, unpinned1, unpinned2, unpinned3, ...
    pinned0, pinned1, pinned2, pinned3)
%CLUSTER Summary of this function goes here
%   Detailed explanation goes here

    index = pins(1,P);

    %Fracture Pins
    addedPins = [];
    [pins, addedPins] = recursive(pins, index, max_rows, max_col, addedPins);
    pins(8,(pins(8,:)==2)) = 1;

    %Find the smallest neighbor of newly broken nodes
    [LNBindex] = smallestNB(addedPins, MLVandBL, max_rows, max_col);

    if isempty(LNBindex)
        return;
    end

    Leng = size(pins,2);
    
    i = ceil(index/max_rows);
    j = mod(index-1,max_col)+1;    
    x = i-xcen;
    y = j-ycen;
    centerDist = sqrt(x^2+y^2);

    %1.4
    %Longer chains at larger radiuses
    threshold = thresholdEQN(centerDist, maxRadius);
        
    %Create a new chain of pins at lowest neighbor
    [pins, ~] = expandClus(pins(:,P), pins, MLVandBL, ...
        LNBindex, threshold, max_rows, max_col);

    %Store all new pins
    for q = Leng+1:size(pins,2)

        index = pins(1,q);
        i = ceil(index/max_rows);
        j = mod(index-1,max_col)+1;

        x = i - xcen;
        y = j - ycen;
        currRad = sqrt(x^2+y^2);

        [pins] =createPin(xcen, ycen,i, j,... 
            pins, maxRadius, currRad, index, q);

        pins(10,q) = 1;
        pins(11,q) = pins(11,P);

        %Find the location of the new pin
        section0 = find(unpinned0 == index);
        section1 = find(unpinned1 == index);
        section2 = find(unpinned2 == index);
        section3 = find(unpinned3 == index);

        %Remove the new pin from unpinned section
        if ~isempty(section0)
            unpinned0(section0) = [];
            pinned0 = [pinned0; index];
        elseif ~isempty(section1)
            unpinned1(section1) = [];
            pinned1 = [pinned1; index];
        elseif ~isempty(section2)
            unpinned2(section2) = [];
            pinned2 = [pinned2; index];
        elseif ~isempty(section3)
            unpinned3(section3) = [];
            pinned3 = [pinned3; index];
        end
    end

    %Replot/store the new pins onto the MLVandBL array
    for q = 1:size(pins,2)
        index = pins(1,q);
        i = ceil(index/max_rows);
        j = mod(index-1,max_col)+1;

        if pins(8,q) > 0 %Check if fractured pin
            MLVandBL(i,j) = 4; %Color to Fracture
        else
            MLVandBL(i,j) = 3; %Color to Pin
        end
    end

end

