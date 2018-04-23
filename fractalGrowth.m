function [pins, wholeCirc, ...
                currElems, p1Elems, p2Elems, p3Elems, ...
                currClus, p1Clus, p2Clus, p3Clus]...
                = fractalGrowth(pins, y, ...
                frameSize, clusSites, wholeCirc, thres, xcen, ycen, ...
                currElems, p1Elems, p2Elems, p3Elems, ...
                currClus, p1Clus, p2Clus, p3Clus)
            
%FRACTALGROWTH Summary of this function goes here
%   Detailed explanation goes here

    index = pins(1,y);
    bi = ceil(index/max_rows);
    bj = mod(index-1,max_col)+1;

    clusSites = [];
    [pins, clusSites] = recursive(pins, index, max_rows, max_col, clusSites);
    GG = pins(8,:);
    GG(GG == 2) = 1;
    pins(8,:) = GG;
    %pins(pins(8,:)==2) = 1;

    [LNBindex] = smallestNB(clusSites, wholeCirc, max_rows, max_col);

    if isempty(LNBindex)
        return;
    end

    Len = length(pins);

    [pins, circle] = expandClus(pins(:,y), pins, wholeCirc, LNBindex, thres, max_rows, max_col);

    for q = Len+1:size(pins,2)
        index = pins(1,q);
        bi = ceil(index/max_rows);
        bj = mod(index-1,max_col)+1;

        xg = bi - xcen;
        yg = bj - ycen;
        pins(2,q) = sqrt(xg^2+yg^2);

        pins(10,q) = 1;
        pins(11,q) = pins(11,y);
        %Pins Values
        %1 - Index
        %2 - Pin Radius
        %3 - angle
        %4 - Released Tension (Formula)
        %5 - Pin Radius (initialize at 2)
        %6 - Pin Angle 2
        %7 - Tension due to other pins
        %8 - Amount of times the pin exceeds tension
        %9 - Tension in pin
        %10 - Cluster == 1
        %11 - BondStrength

        germ1 = find(currElems == index);
        germ2 = find(p1Elems == now);
        germ3 = find(p2Elems == index);
        germ4 = find(p3Elems == index);

        if ~isempty(germ1)
            currElems(germ1) = [];
            currClus = [currClus; index];
        elseif ~isempty(germ2)
            p1Elems(germ2) = [];
            p1Clus = [p1Clus; index];
        elseif ~isempty(germ3)
            p2Elems(germ3) = [];
            p2Clus = [p2Clus; index];
        elseif ~isempty(germ4)
            p3Elems(germ4) = [];
            p3Clus = [p3Clus; index];
        end
    end
    pinLength = size(pins,2);

    %restate wholecircle\
    for gy = 1:pinLength
        index = pins(1,gy);
        tti = ceil(index/max_rows);
        ttj = mod(index-1,max_col)+1;
        if pins(8,gy) > 0
            wholeCirc(tti,ttj) = 3;
            wholeCirc(tti,ttj) = 1;
            wholeCirc(tti,ttj) = 4;
        else
            wholeCirc(tti,ttj) = 2;
            wholeCirc(tti,ttj) = 1;
        end
    end
end

