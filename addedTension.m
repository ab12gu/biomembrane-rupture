function [force] = addedTension(expRadius, pins, maxRadius, pinLength)
%ADDEDTENSION 
    %Calculates the additional tnesion due to pins that are closer
    %to the center of the spread out lipid bilayer

        pinRad = expRadius;
        pinAngle = pins(3,pinLength);
        pinAngle2 = pins(6,pinLength);
        force = 0;

        %Calculate the added tension from closer pins
        for h = 1:pinLength-1
            pinRadz = pins(2,h);
            radDiff = pinRad - pinRadz;

            if pinAngle <= pi/2 && pinAngle > -pi/2
                pinAngz = pins(3,h);
                angDiff = abs(pinAngle - pinAngz);
            else
                pinAngz = pins(6,h);
                angDiff = abs(pinAngle2 - pinAngz);
            end

            if angDiff >= pi/2
                continue;
            end

            %The further the radius (btwn pins), the higher tension
            force = pinForce(radDiff, maxRadius, angDiff, force);
        end

end

