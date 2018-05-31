function [expRadius, new_circle] = Bilayer_Expansion(iteration, center, new_circle)
%EXPANSION 
%   This function expands the current bilayer by a single node radially
    
    expFrame = iteration-1;  %Frame size need for expansion
    [xx, yy] = meshgrid(1:expFrame);
    
    expRadius = iteration/4;  %Radius of Expansion
    origin = iteration/2;
    
    %Fill in all points that are less than radius
    A = sqrt((xx-origin).^2+(yy-origin).^2) <= expRadius;
    
    %High + Low = size of expansion frame
    high = floor(expFrame/2)-1;
    low = ceil(expFrame/2); 
    
    %Resize Circle Frame
    new_circle(center-high:center+low, center-high:center+low) = A;
    %new_circle is the new output

end

