function [pins, tension] = tensionCalc(pins, expRadius, k, pinRad, ...
            maxRadius, P)
%TENSIONCALC Calculates the tension in the pin

    %Calculate the tension due to pins closer to center(7)
    pinToPinTension = pins(7,P)*expRadius^2*k; %k is equal to zero

    %Calculate the total tension in pin
    tension = (expRadius-pinRad)/maxRadius + pinToPinTension;

    %Check if pin is fractured(9) and has replaced its tension with
        %the source fracture's tension
    %Check if pin is not a cluster(10)
    if pins(9,P) > tension && pins(10,P) ~= 1
         %Retain previous tension
         return
    end

    pins(9,P) = tension;
end


%%%%Pin-Row Array Values:
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
