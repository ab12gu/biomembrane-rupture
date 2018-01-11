function[pins, invadedPins, MLVandBL, n] = ...
    circularGrowth(invadedPins, frameSize, tension, pinThres,...
    percReleased, pins, P, xcen, ycen,  MLVandBL, max_rows, max_col)
%CIRCULARGROWTH Summary of this function goes here
%   Detailed explanation goes here

            pinarray = zeros(frameSize);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Release Tension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pins(4,P) = (tension - pinThres)*(1+percReleased/100);
            %pinRadius/currRadius
                        
            %Amount of times pin exceeds tension limit
            pins(8, P) = pins(8, P) + 1;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Expansion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Call Index
            index = pins(1,P);
            i = ceil(index/max_rows);
            j = mod(index-1,max_col)+1;
            
            %Expansion by even numbers
            n = 2*pins(8,P)+pins(5,P) + 2;
            pins(3,P) = n;
            
            %Total size of circle ouput 
            expFrame = n-1;
            [xx, yy] = meshgrid(1:expFrame);
            currRad = n/4;

            %Fill points that are less than radius
            A = sqrt((xx-n/2).^2+(yy-n/2).^2) <= currRad;
            
            %Resize A to correct output ::::: 
            radi = floor((n-1)/2);
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
            
            TempCc = zeros(2,pinGrowth);
            
            %List of expanded nodes: TempCc
            %Index: (1,:) || Radius: (2,:)
            for Rr = 1: pinGrowth
                ni = w(Rr,1);
                nj = w(Rr,2);
                MLVandBL(ni,nj) = 4;
                indexRr = (ni-1)*max_rows+nj;
                check = find(pins(1,:) == indexRr);
                
                if ~isempty(check) && pins(8,check) == 0
                    TempCc(1,Rr) = indexRr;
                    x = ni - xcen;
                    y = nj - ycen;
                    TempCc(2,Rr) = sqrt(x^2+y^2);
                    pins(8,check) = pins(8,P);
                end
            end
            
            tempCirc = [invadedPins TempCc];
            tempCirc( :, all(~tempCirc,1) ) = [];
            invadedPins = tempCirc;

            %Sort Nodes
            [d1, d2] = unique(invadedPins(1,:));
            pinCirc2 = invadedPins(:,d2);
            [d1, d2] = sort(pinCirc2(2,:));
            invadedPins = pinCirc2(:,d2);  
end

