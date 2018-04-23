function [pins, invadedPins, MLVandBL, v] = ...
    growthChain(pins, invadedPins, frameSize, ...
    max_rows, max_col, pinRad, xcen, ycen, MLVandBL, pinProb,... 
    clusProb, maxRadius, expRadius, videoSave, MLVnodes, P, v)
%GROWTHCHAIN Iteratively expands invaded pins
%   Detailed explanation goes here

    while(1)
        if isempty(invadedPins)
            break;
        end
        pinarray = zeros(frameSize);

        index = invadedPins(1,1);
        q = find(index == pins(1,:));
        
        pins(11,q) = pins(11,P);
        pins(8,q) = pins(8,q) + 1;

        i = ceil(index/max_rows);
        j = mod(index-1,max_col)+1;

        for currPinRad = 4:2:pinRad
            %total size of circle ouput 
            expFrame = currPinRad-1;
            [xx, yy] = meshgrid(1:expFrame);

            currRad = currPinRad/4;

            %Fill points that are less than radius
            A = sqrt((xx-currPinRad/2).^2+(yy-currPinRad/2).^2) <= currRad;

            radi = floor((currPinRad-1)/2);
            pinarray(i-radi:...
                i+radi,...
                j-radi:...
                j+radi) = A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Store Expanded Nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            clear w;
            [w(:,1), w(:,2)] = find(pinarray);
            pinGrowth = length(w);

            invadedPins2 = zeros(2,pinGrowth);

            %List of expanded nodes: TempCc
            %Index: (1,:) || Radius: (2,:)
            for q = 1: pinGrowth
                ni = w(q,1);
                nj = w(q,2);
                MLVandBL(ni,nj) = 4;
                index2 = (ni-1)*max_rows+nj;
                check = find(pins(1,:) == index2);

                if ~isempty(check) 
                    if pins(8,check) == 0 && pins(10,check) ~= 1
                        invadedPins2(1,q) = index2;
                        x = ni - xcen;
                        y = nj - ycen;
                        invadedPins2(2,q) = sqrt(x^2+y^2);
                        pins(8,check) = pins(8,P);
                    end
                end
            end

            invadedPinsT = [invadedPins invadedPins2];
            invadedPinsT( :, all(~invadedPinsT,1) ) = [];
            invadedPins = invadedPinsT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visual Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            v = output(MLVandBL, pinProb, clusProb, maxRadius, ...
                frameSize, expRadius, videoSave, pins, MLVnodes, v);

        end

        invadedPins(:,1) = [];

        %Sort Nodes
        [d1, d2] = unique(invadedPins(1,:));
        pinCirc2 = invadedPins(:,d2);
        [d1, d2] = sort(pinCirc2(2,:));
        invadedPins = pinCirc2(:,d2);
    end
end

