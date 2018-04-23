function [index] = indxx(i,j,rowLength)
%Solve for the index (Creates a indexed sequence to quickly reference any
%place/value in the matrix.
index = (i-1)*rowLength+j;

end

