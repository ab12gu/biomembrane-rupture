% -------------------------------------------------------------------------
% Circular Growth Properties (Covering all Failure Modes)
% By: Abhay Gupta
% 08/06/2017
% Rules:
% 1. Assign a source size
% 2. Expand source radius by +2 each iteration
% 3. Give a density of pinning within expansion (probability)
% 4. Give invaders a random number between (1 2)
% 5. Check tension of pin, if tension in pin exceed threshold, 
%    break and spread into a random neighbor
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Clear System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

%opengl hardwarebasic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create a video file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = VideoWriter(['C:\Users\ab18g\DELETE.avi']);
v.FrameRate = 0.0001;
open(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Color Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map = [...
    0.00 0.00 0.00;     %Background     (Black)
    %0.85 0.84 0.87;  %Background || 0(1) < x=< 2   (Light Grey)
    1 0.1 0.1;
    %0.30 0.50 0.77;     %Source     || 2    < x=< 3   (Turquoise)
    %0.20 0.30 0.40;     %Bi-layer   || 3    < x=< 4   (Blue)
    0.85 0.0 0.0;     %Invader    || 3 < x=< 4   (Brush Red)
    0.5 0.50 1;     %Pins       || 4 < x=< 5   (Hot Pink)
    %0.90 0.00 0.00;     %Pins       || 4    < x=< 5   (Black)
    0.5 0.0 0.0];            %Fracture   || 5    < x=< inf (White)
%The map increments by 1.

%Colors to match the experimental output:
%0.85 0.10 0.20;     %Invader    || 3 < x=< 4   (Brush Red)
%0.95 0.00 0.34;     %Pins       || 4 < x=< 5   (Hot Pink)

colormap(map)  
%NOTE: output = Array+1
        
%Legend:
%Background         = 0 (or 1)
%Source             = 2
%Displaced Bi-layer = 3
%Pins               = 4
%Fracture           = 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Changable Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make invader ceil for output
maxRadius = 200;

%Density/Probability of pins
prob = 0.001;
factor = 1;
probMin = prob*factor;
factor2 = 1/9
radMin = maxRadius*factor2

%Likliness of pin becoming a cluster
pl = 1;

%Threshold Value for cluster size
threshold = 0.7; %0.35+1;

%Tension Limit for circle expansion
tensionThreshold = 0.6; %0.7;

%Percent Tension Released
released = 8;

%Biasing toward expansion
bias = 1;
decrement = 1;

%Factor of pinToPinTension
k = 0;

%Cluster Tension Parameters:
%Constant
kc = 1;
%Exponential
pc = 1;

%Pin Bond Strength
B1 = 0.6;

%Cluster Bond Added Strength
B2 = 0.6;

%Radius of source
initialRadius  = 6;

%Initial Output Radius 
OutputRadiusInit = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial Calculated Paramters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Maximum area of displaced bi-layer
maxArea = pi*maxRadius^2;

%Frame/Array size
frameSize = maxRadius*4;
max_rows = frameSize;
max_col = frameSize;
PinAmount = floor(prob*maxArea);

%Location of MLV 
xcen = (max_col+1)/2;
ycen = (max_rows+1)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Pre-allocating Memory
CircCurr = zeros(frameSize);
circPrev = zeros(frameSize);
Circum   = zeros(frameSize);
%pins = zeros(8,floor(prob*maxArea));

% Matrix of neighbors
NB = zeros(9,1);

% Amount of pins & invaders
pinLength = 0;
invLength = 0;
pins = [];

%Lowest Neighbor Value for internal percolation
LNB = threshold;

%Number of pin clusters
clusLength = 0;

%Elems is unpinned Lipids
%Clus is pinned Lipids
expPrev     = 0;
currElems   = 0;
currClusLng = [];
p3Elems     = [];
p2Elems     = [];
p1Elems     = [];
p3Clus      = [];
p2Clus      = [];
p1Clus      = [];
currClus    = [];
pinCirc     = [];

%Set the initial MLV radius
%NOTE: 'm' is a 4 times effective radius
start = (initialRadius)*4+1;
minM = start+4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Begin lipid bi-layer spread across Si02 surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m = start:2:frameSize
    
    %total size of circle ouput 
    circleArray = m-1;
    [rr, cc] = meshgrid(1:circleArray);
    
    currRadius = m/4;
    Radius = round(currRadius);
    
    %Fill points that are less than radius
    A = sqrt((rr-m/2).^2+(cc-m/2).^2) <= currRadius;

    %Resize Circle Frame
    CircCurr(floor(max_rows/2)-floor(circleArray/2)+1:...
        floor(max_rows/2)+ceil(circleArray/2),...
        floor(max_col/2)-floor(circleArray/2)+1:...
        floor(max_col/2)+ceil(circleArray/2)) = A;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if m == start
        %Change color of initial source
        %Keep source value at 1
        wholeCirc = CircCurr;
        circPrev = wholeCirc;
        
        %Store invader into invader array
        [f(:,1), f(:,2), f(:,3)] = find(wholeCirc);
        growthSize = length(f);
        
        for h = 1:growthSize
            i = f(h, 1);
            j = f(h, 2);
            index = (i-1)*max_rows+j;
            invLength = invLength + 1;
            C(1,invLength) = index;
        end
        
        clear f;
        continue;
    end
    
    %Store the new nodes and reset previous circle
    Circum = CircCurr - circPrev;
    circPrev = CircCurr;
    
    %Find the Index and assign a value to each new element
    [f(:,1), f(:,2), f(:,3)] = find(Circum);
    growthSize = length(f);
    
    Circum = zeros(frameSize);
        
    %Store the last four expansion's invader elements (not pins)   
    p3Elems = p2Elems;
    p2Elems = p1Elems;
    p1Elems = currElems;
    clear currElems;
    
    %Store the last four expansion's cluster elements (pins)
    p3Clus = p2Clus;
    p2Clus = p1Clus;
    p1Clus = currClus;
    currClus = [];
    
    %set current invader length
    invLength = 0;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find all pinning sites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for h = 1:growthSize
        
        %Find index
        i = f(h, 1);
        j = f(h, 2);
        index = (i-1)*max_rows+j;
                
        %Set to pin density
        if rand(1) < (prob-probMin)/(maxRadius-radMin)*(Radius)+probMin && Radius>radMin %rob %Turn node into a pin          
            %Note pin in new matrix f
                %f3 remains a 1
            
            %Count the additional pins
            pinLength = pinLength + 1;
            
            %Store pin into output
            %Start pin count at 3.
            Circum(i,j) = pinLength+2;

            %Find the angle of the pin
            x = j - xcen;
            y = i - ycen;
            %Angles: -pi -> 0 || 0 -> pi
            pins(3,pinLength) = atan2(y,x);
            %Angles: 0 -> 2pi
            pins(6,pinLength) = mod(atan2(y,x)+2*pi,2*pi);
            
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
            
            %Store the index & pin radius
            pins(1,pinLength) = index;
            pins(2, pinLength) = currRadius;
            
            %Tension released
            pins(4, pinLength) = 0;
                   
            %current radius of pin
            pins(5, pinLength) = 2;
            
            pinRad = currRadius;
            pinAngle = pins(3,pinLength);
            pinAngle2 = pins(6,pinLength);
            force = 0;
            
            %Calculate the added tension from closer pins
            for z = 1:pinLength-1
                pinRadz = pins(2,z);
                radDiff = pinRad - pinRadz;

                if pinAngle <= pi/2 && pinAngle > -pi/2
                    pinAngz = pins(3,z);
                    angDiff = abs(pinAngle - pinAngz);
                else
                    pinAngz = pins(6,z);
                    angDiff = abs(pinAngle2 - pinAngz);
                end
                
                if angDiff >= pi/2
                    continue;
                end
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Effect of pins on pins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %The further the radius (btwn pins), the higher tension
                force = pinForce(radDiff, maxRadius, angDiff, force);
            end
            
            %disp(force);
            pins(7, pinLength) = force;
            
            %Amount of times pin exceeds tension limit
            pins(8, pinLength) = 0;
            
            %Pin Tension
            pins(9, pinLength) = 0;

            %Check if site becomes a cluster
            if rand(1) < pl
                currClus = [currClus;index];
                pins(10,pinLength) = 1;
                pins(11,pinLength) = B1+B2;
            else
                pins(10,pinLength) = 0;
                pins(11,pinLength) = B1;
            end
            
            %Label
        %Else, add a regular invader node
        else
            %if not a pin, change name to not pin and add rand number
            Circum(i,j) = rand*bias+1;
            
            %store the index of the invader
            invLength = invLength + 1;
            currElems(invLength, 1) = index;
            
        end
    end   
        
    %Total invader nodes from 4 iterations
    pcElems = [p3Elems; p2Elems; p1Elems; currElems];
    pcClus = [p3Clus; p2Clus; p1Clus; currClus];    
    rh = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spread clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    while(rh ~= length(pcClus)) 
                
        if m <= minM
            break;
        end
        
        %Increment
        rh = rh + 1;

        %index of cluster
        indexClus = pcClus(rh);

        i = ceil(indexClus/max_rows);
        j = mod(indexClus-1,max_col)+1;
        
        %check neighbors if they are lower than threshold
        %if so, check if within pcElems
        %if so, expand and store as new pin and cnew cluster
        %increment an overall cluster addition
        xi = i-xcen;
        xj = j-ycen;
        CR = sqrt(xi^2+xj^2);
        
        %1.4
        %Longer chains at larger radiuses
        thres = 0.70*log2(0.5*(CR/maxRadius)+1)+1;
        %2.3*log10(0.5*(CR/maxRadius)+1)+1;
        
        for mv = 1:3
            for nv = 1:3
                
                %jerry = wholeCirc;
                %find the neighbor's value
                nbi = i+(mv-2);
                nbj = j+(nv-2);
                indxnb = (nbi-1)*max_rows+nbj;
                
                %Check if element is invader for pc and curr
                indxpc = find(pcElems == indxnb);
                indxcurr = find(currElems == indxnb);
                
                if  isempty(indxcurr)
                    T = wholeCirc(nbi, nbj);
                else
                    T = Circum(nbi ,nbj);
                end

                %if neighbor value is a pin value, store as simple invader
                if T >= 2 || T <= 1
                    T = 2;
                end
                                
                %Check if the neighbor is in the element range, create pin
                %and below threshold
                if T < 2 && ~(isempty(indxpc)) && T < thres
                    
                    indxp1E  = find(p1Elems == indxnb);
                    indxp2E  = find(p2Elems == indxnb);
                    indxp3E  = find(p3Elems == indxnb);
                    %Add pin to circle
                    
                    %Count the additional pins
                    pinLength = pinLength + 1;
            
                    %Store pin into output
                    %Start pin count at 3.
                    %Delete invader Element, add cluster Element
                    if ~isempty(indxcurr)
                        Circum(nbi,nbj) = pinLength+2;
                        currElems(indxcurr) = [];
                        currClus = [currClus; indxnb];
                    else
                        wholeCirc(nbi, nbj) = pinLength +2;
                        if ~isempty(indxp1E)
                            p1Elems(indxp1E) = [];
                            p1Clus = [p1Clus; indxnb];
                        elseif ~isempty(indxp2E)
                            p2Elems(indxp2E) = [];
                            p2Clus = [p2Clus; indxnb];
                        elseif ~isempty(indxp3E)
                            p3Elems(indxp3E) = [];
                            p3Clus = [p3Clus; indxnb];
                        end
                    end
                    
                    %add new cluster to cluster list
                    pcClus = [pcClus; indxnb];
                    
                    %remove node from both curr & pc elem
                    pcElems(indxpc) = [];
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pin Storage & Tension Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %Store pin in pin array
                    pins(1,pinLength) = indxnb;
                    index = indxnb;
                    i = ceil(index/max_rows);
                    j = mod(index-1,max_col)+1;
                                        
                    %Find the angle of the pin
                    %The distance from pin to center:
                    x = i - xcen;
                    y = j - ycen;
                    currRad = sqrt(x^2+y^2);
                    pinRad = currRad;

                    %Angles: -pi -> 0 || 0 -> pi
                    pinAngle = atan2(y,x);
                    pins(3,pinLength)= pinAngle;
                    %Angles: 0 -> 2pi
                    pinAngle2 = mod(atan2(y,x)+2*pi,2*pi);
                    pins(6,pinLength) = pinAngle2;
                    
                    %Store the index & pin radius
                    pins(1,pinLength) = index;
                    pins(2, pinLength) = currRad;
            
                    %Tension released
                    pins(4, pinLength) = 0;
                    %Current radius of pin
                    pins(5, pinLength) = 2;
                    %Initialize force from other pins
                    force = 0;
            
                    %Calculate the added tension from closer pins
                    for z = 1:pinLength-1
                        pinRadz = pins(2,z);
                        radDiff = pinRad - pinRadz;

                        if pinAngle <= pi/2 && pinAngle > -pi/2
                            pinAngz = pins(3,z);
                            angDiff = abs(pinAngle - pinAngz);
                        else
                            pinAngz = pins(6,z);
                            angDiff = abs(pinAngle2 - pinAngz);
                        end

                        if angDiff >= pi/2
                            continue;
                        end
                        %The further the radius (btwn pins), the higher tension
                        force = pinForce(radDiff, maxRadius, angDiff, force);
                    end

                    %Tension on pin (disp force)
                    pins(7, pinLength) = force;

                    %Amount of times pin exceeds tension limit
                    Loc = find(indexClus == pins(1,:));
                    pins(8, pinLength) = pins(8,Loc);
                    
                    %Amount of tension
                    pins(9, pinLength) = 0;
                    
                    %Store as a cluster
                    pins(10,pinLength) = 1;
                    
                    Loc = find(indexClus == pins(1,:));
                    Stretch = 7/8;
                    
                    %Store the new bond strength
                    pins(11,pinLength) = pins(11,Loc)*Stretch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visual of expansion in pin clusters (jerry)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
%                     jerry = jerry + Circum; 
%                     
%                     for pn = 1:size(pins,2)
%                        if pins(8,pn) == 0
%                             indxpin = pins(1,pn);
%                             ki = ceil(indxpin/max_rows);
%                             kj = mod(indxpin-1,max_col)+1;
%                             
%                             jerry(ki,kj) = 2.7;
%                        end
%                     end
%                     
%                     image(ceil(jerry)+1)
%                     axis image
%                     h = colorbar;  
%                     set(h, 'XTick', [1.5, 2.5, 3.5, 4.5, 5.5])
%                     set(h,'XTickLabel', {   'Background',...
%                                             'Source (MLV)',...
%                                             'Displaced \n /n BiLayer',...
%                                             'Pinning Sites',...
%                                             'Fractured Area'})
%                     title(h, 'Legend')                                  
%                     %frame=getframe(gcf);
                end
            end
        end
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add new nodes to output (wholeCirc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bias = bias*decrement;
    wholeCirc = wholeCirc + Circum; 
    
    for pn = 1:size(pins,2)
        indxpin = pins(1,pn);
        ki = ceil(indxpin/max_rows);
        kj = mod(indxpin-1,max_col)+1;
        
        if pins(8,pn) == 0
            wholeCirc(ki,kj) = 3;
        else
            wholeCirc(ki,kj) = 4;
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Recalculate Tension (after each iteration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pinCirc = [];
    tempCirc = [];

    y = 0;
    for y = 1:pinLength
            
        pinRad = pins(2,y);
        pinAngle = pins(3,y);
        pinAngle2 = pins(6,y);
                                                                                 
        %Note: everything above pin value is has a smaller or equal radius
     
        %Calculate Added Tension in pin due to closer pins
        %Only needs to be recalculated if intricate equation
        %May need to be added now to incorporate change in tension due to 
        %fracture of the inner pins
        
        %for z = 1:pinLength-1
        %    pinRadz = pins(2,z);
        %    radDiff = pinRad - pinRadz;
        %    
        %    if pinAngle <= pi/2 && pinAngle > -pi/2
        %        pinAngz = pins(3,z);
        %        angDiff = abs(pinAngle - pinAngz);
        %    else
        %        pinAngz = pins(6,z);
        %        angDiff = abs(pinAngle2 - pinAngz);
        %    end
        %end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tension due to other pins (dynamic with radius) & Total tension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %Emily = pins(10,y); 
        pinToPinTension = pinToPinTens(pins(7,y),currRadius,k);
        pinRR = pins(2,y);
        %clusterTension = Emily*kc*(1-((currRadius-pinRR)/maxRadius)^pc);
        tensionReleased = pins(4,y);
        %Total Tension
        tension = (currRadius- pinRad)/maxRadius - tensionReleased...
            + pinToPinTension; %+ clusterTension - Emily/2;
        if pins(9,y) > tension && pins(10,y) ~= 1
            tension = pins(9,y);
        end
        pins(9,y) = tension;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fractal Growth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if tension > pins(11,y) && pins(10,y) > 0
            if pins(8,y) == 0
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
                    break;
                end
                    
                Len = length(pins);
                
                [pins, circle] = expandClus(pins(:,y), pins, wholeCirc, LNBindex, thres, max_rows, max_col);
                
                for q = Len+1:length(pins)
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
        end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CircularGrowth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Check Tension & if not cluster
        if tension > pins(11,y) && pins(10,y) == 0
            %Initialize Board
            pinarray = zeros(frameSize);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Release Tension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pins(4,y) = (tension - tensionThreshold)*(1+released);
            %pinRadius/currRadius
                        
            %Amount of times pin exceeds tension limit
            pins(8, y) = pins(8, y) + 1;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Expansion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Call Index
            index = pins(1,y);
            bi = ceil(index/max_rows);
            bj = mod(index-1,max_col)+1;
            
            %Expansion by even numbers
            n = 2*pins(8,y)+pins(5,y) + 2;
            pins(3,y) = n;
            
            %Total size of circle ouput 
            circleArray = n-1;
            [rr, cc] = meshgrid(1:circleArray);
            currRad = n/4;

            %Fill points that are less than radius
            A = sqrt((rr-n/2).^2+(cc-n/2).^2) <= currRad;
            
            %Resize A to correct output ::::: 
            radi = floor((n-1)/2);
            pinarray(bi-radi:...
                bi+radi,...
                bj-radi:...
                bj+radi) = A;
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Store Expanded Nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            clear w;
            [w(:,1), w(:,2)] = find(pinarray);
            pinGrowth = length(w);
            expSize = size(pinCirc,2);
            
            TempCc = zeros(2,pinGrowth);
            
            %List of expanded nodes: TempCc
            %Index: (1,:) || Radius: (2,:)
            for Rr = 1: pinGrowth
                ni = w(Rr,1);
                nj = w(Rr,2);
                wholeCirc(ni,nj) = 4;
                indexRr = (ni-1)*max_rows+nj;
                check = find(pins(1,:) == indexRr);
                
                if ~isempty(check) && pins(8,check) == 0
                    TempCc(1,Rr) = indexRr;
                    xi = ni - xcen;
                    xj = nj - ycen;
                    TempCc(2,Rr) = sqrt(xi^2+xj^2);
                    pins(8,check) = pins(8,y);
                end
            end
            
            tempCirc = [pinCirc TempCc];
            tempCirc( :, all(~tempCirc,1) ) = [];
            pinCirc = tempCirc;

            %Sort Nodes
            [d1, d2] = unique(pinCirc(1,:));
            pinCirc2 = pinCirc(:,d2);
            [d1, d2] = sort(pinCirc2(2,:));
            pinCirc = pinCirc2(:,d2);         
        end
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Expand the invaded pins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        if ~isempty(pinCirc)
            while(1)
                if isempty(pinCirc)
                    break;
                end
                pinarray = zeros(frameSize);
                
                index = pinCirc(1,1);
                xx = find(pins(1,:) == indexRr);
                
                pins(4,xx) = (tension - tensionThreshold)*(1+released);
                pins(8,xx) = pins(8,xx) + 1;

                bi = ceil(index/max_rows);
                bj = mod(index-1,max_col)+1;

                for ncurr = 4:2:n
                    %total size of circle ouput 
                    circleArray = ncurr-1;
                    [rr, cc] = meshgrid(1:circleArray);

                    currRad = ncurr/4;

                    %Fill points that are less than radius
                    A = sqrt((rr-ncurr/2).^2+(cc-ncurr/2).^2) <= currRad;

                    radi = floor((ncurr-1)/2);
                    pinarray(bi-radi:...
                        bi+radi,...
                        bj-radi:...
                        bj+radi) = A;
                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Store Expanded Nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
                    clear w;
                    [w(:,1), w(:,2)] = find(pinarray);
                    pinGrowth = length(w);
                    expSize = size(pinCirc,2);

                    TempCc = zeros(2,pinGrowth);

                    %List of expanded nodes: TempCc
                    %Index: (1,:) || Radius: (2,:)
                    for Rr = 1: pinGrowth
                        ni = w(Rr,1);
                        nj = w(Rr,2);
                        wholeCirc(ni,nj) = 4;
                        indexRr = (ni-1)*max_rows+nj;
                        check = find(pins(1,:) == indexRr);

                        if ~isempty(check) && pins(8,check) == 0 && pins(10,check) ~= 1
                            TempCc(1,Rr) = indexRr;
                            xi = ni - xcen;
                            xj = nj - ycen;
                            TempCc(2,Rr) = sqrt(xi^2+xj^2);
                            pins(8,check) = pins(8,y);
                        end
                    end
                    
                    tempCirc = [pinCirc TempCc];
                    tempCirc( :, all(~tempCirc,1) ) = [];
                    pinCirc = tempCirc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Show All Pins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for node = 1:pinLength
                        index = pins(1,node);
                        i = ceil(index/max_rows);
                        j = mod(index-1,max_col)+1;
                        if pins(8,node) > 0
                            if pins(10,node) ~= 1
                                wholeCirc(i,j) = 3;
                            else
                                wholeCirc(i,j) = 4;
                            end
                        else
                            wholeCirc(i,j) = 3;
                        end
                    end

                    for node = 1:length(C)
                        index = C(1,node);
                        i = ceil(index/max_rows);
                        j = mod(index-1,max_col)+1;
                        wholeCirc(i,j) = 1;
                    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visual Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    image(ceil(wholeCirc)+1)
                    axis image

                    title({['Pin Density: ',  num2str(prob*100), '%'];...
                        ['Cluster Density: ', num2str(pl*100), '%' ]});
                    set(gca,'xtick',[]);
                    set(gca,'YTick',[]);
                    xlabel({['Maxium Radius:' num2str(maxRadius), ' || Frame Size:' ...
                        num2str(frameSize), ':', num2str(frameSize)];...
                        [' Radius of Displaced Bi-layer: ', ...
                        num2str(Radius)]});

                    %Frame Settings
                    width= 850;
                    height= 600;
                    set(gcf,'position',[300,50,width,height])

                    %Colorbar Settings:
                    h = colorbar;  
                     set(h, 'XTick', [1.5, 2.5, 3.3, 3.7, 4.5, 5.5])
                     set(h,'XTickLabel', {   'Background',...
                                             'Source (MLV)',...
                                             'Bi-Layer',...
                                             'Displaced',...
                                             'Pinning Sites',...
                                             'Fractured Area'})
                    title(h, '                    Colormap Legend')

                    frame=getframe(gcf);
                    writeVideo(v,frame);

                end
                
                pinCirc(:,1) = [];

                %Sort Nodes
                [d1, d2] = unique(pinCirc(1,:));
                pinCirc2 = pinCirc(:,d2);
                [d1, d2] = sort(pinCirc2(2,:));
                pinCirc = pinCirc2(:,d2);
            end
        end
    end
    
    pinLength = size(pins,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Show All Pins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for node = 1:pinLength
        index = pins(1,node);
        i = ceil(index/max_rows);
        j = mod(index-1,max_col)+1;
        if pins(8,node) > 0
            if pins(10,node) ~= 1
                wholeCirc(i,j) = 3;
            else
                wholeCirc(i,j) = 4;
            end
        else
            wholeCirc(i,j) = 3;
        end
    end

    for node = 1:length(C)
        index = C(1,node);
        i = ceil(index/max_rows);
        j = mod(index-1,max_col)+1;
        wholeCirc(i,j) = 1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visual Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Radius > OutputRadiusInit
        image(ceil(wholeCirc)+1)
        axis image

        title({['Pin Density: ',  num2str(prob*100), '%'];...
            ['Cluster Density: ', num2str(pl*100), '%' ]});
        set(gca,'xtick',[]);
        set(gca,'YTick',[])
        xlabel({['Maxium Radius:' num2str(maxRadius), ' || Frame Size:' ...
            num2str(frameSize), ':', num2str(frameSize)];...
            [' Radius of Displaced Bi-layer: ', ...
            num2str(Radius)]});

        %Frame Settings
        width= 850;
        height= 600;
        set(gcf,'position',[300,50,width,height])

        %Colorbar Settings:
        h = colorbar;  
         set(h, 'XTick', [1.5, 2.5, 3.3, 3.7, 4.5, 5.5])
         set(h,'XTickLabel', {   'Background',...
                                 'Source (MLV)',...
                                 'Bi-Layer',...
                                 'Displaced',...
                                 'Pinning Sites',...
                                 'Fractured Area'})
        title(h, '                    Colormap Legend')

        frame=getframe(gcf);
        writeVideo(v,frame);
    end
    
    clear f
end

close(v);
