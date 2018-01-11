% -------------------------------------------------------------------------
% Simulation of Cell Fracture
% By: Abhay Gupta
% 08/06/2017
% Rules:
% 1. Begin with an initial MLV size
% 2. Slowly spread the lipid-bilayer radially
% 3. Strategically/Probabilistically assign pins across expansion
% 4. Strategically/Probabilistically assign cluster of pins across
% expansion
% 5. Define a bond strength to each pin
% 6. Break bonds if the tension in pin exceeds bond strength of pin
% 7. End expansion once the diamater of expansion is half the frame size
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a Video file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
videoSave = 0;

if videoSave == 1
    v = VideoWriter([...
       'C:\Users\ab18g\Google Drive\Research\Spatial Simulation\Video_Output\Theory2\' ...
       'FrameRateCheck2.avi']);
    v.FrameRate = 10; %FrameRate = Frames/Second
    open(v);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Color Settings 
%%% Range from 0 -> 1
%%% [Red Green Blue] Variable        Range          Color
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmap = [...
    0.00 0.00 0.00;  %Background   || -inf <x=< 2   (Black)
    1.00 0.10 0.10;  %MLV          || 2    <x=< 3   (Light Red)
    0.85 0.00 0.00;  %Bilayer      || 3    <x=< 4   (Brush Red)
    0.90 0.00 0.00;  %Pinned BL    || 4    <x=< 5   (Red)
    0.50 0.00 0.00]; %Fractured BL || 5    <x=< inf (Dark Red)
colormap(cmap)  


%NOTE: The plotted value is Array Value + 1

%More Colors:
    %0.95 0.00 0.34;     (Hot Pink)
    %0.30 0.50 0.77;     (Turquoise)
    %0.20 0.30 0.40;     (Blue)
    
%Legend for code:
    %Background = 0 (or 1)
    %MLV        = 1
    %Bilayer    = 2
    %Pinned BL  = 3
    %Fracture   = 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Paramaterized Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxRadius = 200; %Maximum radius for bilayer growth
pinProb = 0.001; %Probability of expanded bilayer nodes turning into pins

%%Nonlinear probability across BL (increase pins probability further
    %%away from center
    factor = 1;
    minProb = pinProb*factor;
    factor2 = 1/9;
    radMin = maxRadius*factor2;

clusProb = 0.8;         %Probability of new pin sites turning into cluster of pins
clusThres = 0.7;        %Threshold Value for cluster expansion
pinThres = 0.6;         %Tension Limit for single pin expansion
percReleased = 8;       %Percent Tension Released after node breaks
singleBond = 0.6;       %Pin Bond Strength
doubleBond = 1.2;       %Additional Cluster Bond Strength (Double Bonds)
stretch = 7/8;          %The amount of stretch each new cluster bond has || Lower bond strength
MLVradius  = 6;         %MLV Radius
minOutputRadius = 100;  %The minimum radius required for image to be graphed
k = 0;                  %Percentage of how much other closer pins affect further pins

%Cluster Tension Parameters:
kc = 1; %Constant
pc = 1; %Exponential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Calculated Paramters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxAreaBL = pi*maxRadius^2; %Maximum Area of the Bilayer 

%Frame/Array size
frameSize = maxRadius*4;
center = frameSize/2;
max_rows = frameSize;
max_col = frameSize;

%Center of MLV 
xcen = (max_col+1)/2;
ycen = (max_rows+1)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Pre-allocating Memory
circCurr = zeros(frameSize);
circPrev = zeros(frameSize);
expansion   = zeros(frameSize);

% Amount of pins & invaders
pinLength = 0;
pins = [];

%Lowest Neighbor Value for internal percolation
LNB = clusThres;

%Expanded unpinned BL nodes for the last 3 expansions
unpinned3    = [];
unpinned2    = [];
unpinned1    = [];
unpinned0    = [];

%Expanded pinned BL nodes for the last 3 expansions
pinned3      = [];
pinned2      = [];
pinned1      = [];
pinned0      = [];

%Set the initial MLV paramters
%NOTE: 'm' is a 4 times the current radius of the bilayer
start = (MLVradius)*4+1;
minIteration = start+4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin lipid bi-layer spread across Si02 surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iteration = start:2:frameSize    
    
    %Recomputes array for the expansion of the bilayer
    [expRadius, circCurr] = expansionBL(iteration, center, circCurr); 
     
    %Creates the initial MLV array
    if iteration == start
        [circCurr, MLVandBL, MLVnodes] = MLV(circCurr, max_rows);         
        continue;
    end
    
    %Stores the new circumference of nodes into expansion
    expansion = circCurr - circPrev; %
    circPrev = circCurr; %Reset Previous Circle to Current Circle
    
    %Find the the new indexes for each node
    [i, j] = find(expansion);
    growthSize = length(i);
    
    expansion = zeros(frameSize);
        
    %Store the last four expansion's unpinned elements   
    unpinned3 = unpinned2;
    unpinned2 = unpinned1;
    unpinned1 = unpinned0;
    unpinned0 = [];
    
    %Store the last four expansion's pinned elements
    pinned3 = pinned2;
    pinned2 = pinned1;
    pinned1 = pinned0;
    pinned0 = [];
    
    %set current invader length
    invLength = 0;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find all pinning sites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for h = 1:growthSize
        
        index = indxx(i(h), j(h), max_rows); 
        %A Non-linear probability of pinning occuring radially
            %expRadius is the changing input variable
        pinDensityValue = (pinProb-minProb)/(maxRadius-radMin)...
            *(expRadius)+minProb; %expRadius: current radius of the BL
                        
        %Save each new pin site into the pin array (pins)
            %radMin is the minimum radius needed for pin sites to form
        if rand(1) < pinDensityValue && expRadius > radMin 

            %Store pin into output array
            expansion(i(h),j(h)) = 3; %+2 for output value
          
            pinLength = size(pins,2)+1;
            
            [pins] = createPin(xcen, ycen, i(h),... 
                j(h), pins, maxRadius, expRadius, index, pinLength);
             
            %Check if site becomes a cluster type pin
            if rand(1) < clusProb
                pinned0 = [pinned0;index]; %Cluster
                pins(10,pinLength) = 1;
                pins(11,pinLength) = doubleBond;
            else
                pins(10,pinLength) = 0; %Not cluster, single pin
                pins(11,pinLength) = singleBond;
            end

        %Else, store as regular expanded BL node
        else
            %if not a pin, change name to not pin and add rand number
            expansion(i(h),j(h)) = rand+1;
            unpinned0 = [unpinned0; index]; %Store as new BL 
        end
    end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spread pinning chains (clusters of pins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Total invader nodes from 4 iterations
    unpinnedNodes = [unpinned3; unpinned2; unpinned1; unpinned0];
    pinnedNodes = [pinned3; pinned2; pinned1; pinned0];    
    pinnedLength = 0;
    clear i j;
                          
    [pins, expansion, MLVandBL, threshold, ...
        pinned0, pinned1, pinned2, pinned3, pinnedNodes, ...
        unpinned0, unpinned1, unpinned2, unpinned3, unpinnedNodes] = ...
        clusterSpread(iteration, minIteration, xcen, ycen, pinnedNodes, unpinnedNodes, ...
        maxRadius, MLVandBL, expansion, unpinned0, unpinned1, ...
        unpinned2, unpinned3, pinned0, pinned1, pinned2, pinned3, stretch, ...
        pinnedLength, pins, max_rows, max_col);
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add new nodes to output (MLVandBL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MLVandBL = MLVandBL + expansion; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recalculate Tension (after each iteration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    invadedPins = [];
    
    for P = 1:pinLength
        pinRad = pins(2,P);
        [pins, tension] = tensionCalc(pins, expRadius, k, pinRad, ...
            maxRadius, P);
                            
        bondStrength = pins(11,P);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fractal Fracture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %If the tension exceeds bond strength 
         %And is a chained pin (cluster)(10)
          %And is not previously broken: (8)
            %Do Fractal Fracture       
        if tension > bondStrength && pins(10,P) > 0 && ~(pins(8,P))
                
            [pins, LNBindex, MLVandBL] = ...
            cluster(pins, P, max_rows, max_col, MLVandBL, ...
                threshold, ycen, xcen, maxRadius, ...
                unpinned0, unpinned1, unpinned2, unpinned3, ...
                pinned0, pinned1, pinned2, pinned3);
        
            if isempty(LNBindex)
                continue;
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CircularGrowth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Check Tension & if not cluster
        if tension > pins(11,P) && pins(10,P) == 0
            
            [pins, invadedPins, MLVandBL, n] = ... 
                circularGrowth(invadedPins, frameSize, tension, pinThres, ...
                percReleased, pins, P, xcen, ycen,  MLVandBL, ...
                max_rows, max_col); 
        end
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Expand the invaded pins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Check if any pins circular fractured into previously pinned
        %locations
        if ~isempty(invadedPins)
            [pins, invadedPins, MLVandBL] = growthChain(pins, ...
                invadedPins, frameSize, tension, pinThres, percReleased, ...
                max_rows, max_col, n, xcen, ycen, MLVandBL, pinProb, ...
                clusProb, maxRadius, expRadius, videoSave, MLVnodes, P);
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visual Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if expRadius > minOutputRadius
        output(MLVandBL, pinProb, clusProb, maxRadius, frameSize, ...
            expRadius, videoSave, pins, MLVnodes);
    end
    
    clear f
end

if videoSave == 1
    close(v);
end

%%END of File