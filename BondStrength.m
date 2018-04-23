function [bondStrength] = BondStrength(Length,stretch,initialBondStr)
%BONDSTRENGTH Determines the bond strength of the pin

    bondStrength = (stretch)^(Length)*initialBondStr*(Length+1);

end

