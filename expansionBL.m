function [expRadius, CircCurr] = expansionBL(miniframe, center, CircCurr)
%EXPANSION 
%   This function expands the current bilayer by a single node radially
    %miniframe = iteration
    
    expFrame = miniframe-1;  %Frame size need for expansion
    [xx, yy] = meshgrid(1:expFrame);
    
    expRadius = miniframe/4;  %Radius of Expansion
    origin = miniframe/2;
    
    %Fill in all points that are less than radius
    A = sqrt((xx-origin).^2+(yy-origin).^2) <= expRadius;
    
    %High + Low = size of expansion frame
    high = floor(expFrame/2)-1;
    low = ceil(expFrame/2); 
    
    %Resize Circle Frame
    CircCurr(center-high:center+low, center-high:center+low) = A;
    %CircCurr is the new output

end

