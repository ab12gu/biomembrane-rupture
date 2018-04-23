function [pineven, pinodd] = nodeAngle(i, j, xcen, ycen)
%PINANGLE Determines the angle of the pin from the center
    %Stores the angle in within 2 seperate periods
        x = j - xcen;
        y = i - ycen;
        
        %Angles: -pi -> 0 || 0 -> pi
        pineven = atan2(y,x);
        
        %Angles: 0 -> 2pi
        pinodd = mod(atan2(y,x)+2*pi,2*pi);

end

