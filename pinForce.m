function [output] = pinForce(radDiff, maxRadius, angDiff, force)
%Function Name: pinForce
%   Function that defines the force caused by pins near center

output = (1-radDiff/maxRadius)*cos(angDiff)+force;

end

