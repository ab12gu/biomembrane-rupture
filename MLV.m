function [CircCurr, wholeCirc, MLVnodes] = MLV(CircCurr, max_rows)
%SOURCE This function finds the index values of the initial MLV

    %Change color of initial source
    %Keep source value at 1
    wholeCirc = CircCurr;

    %Store BL expansion into BL expansion array
    [i, j] = find(CircCurr);
    sourceSize = length(i);

    for h = 1:sourceSize
        index = indxx(i(h), j(h), max_rows);
        MLVnodes(1,h) = index;
    end
end

