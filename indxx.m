function [index] = indxx(i,j,max_row_length)
%Solve for the index (Creates a indexed sequence to quickly reference any
%place/value in the matrix)
index = (i-1)*max_row_length+j;

end

