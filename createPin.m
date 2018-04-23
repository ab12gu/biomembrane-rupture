function [pins] = ...
    createPin(xcen, ycen, i, j, pins, maxRadius, ...
        expRadius, index, pinLength)
%CREATEPIN Stores pin into an array

    %Find the angle of the pin
    [pineven, pinodd] = nodeAngle(i, j, xcen, ycen);

    pins(1,pinLength) = index;
    pins(2,pinLength) = expRadius;
    pins(3,pinLength) = pineven; % -pi < angle <  pi
    pins(4,pinLength) = 0;
    pins(5,pinLength) = 2;
    pins(6,pinLength) = pinodd;  % 0   < angle < 2pi

    %Calculate the added tension from pins closer to the center of
    %the spreaded Lipid Bilayer

    [force] = addedTension(expRadius, pins, maxRadius, pinLength);
    pins(7, pinLength) = force;
    pins(8, pinLength) = 0;
    pins(9, pinLength) = 0;
                
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