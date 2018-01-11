function [LNBindex] = smallestNB(addedPins, MLVandBL, max_rows, max_col)
%Find the smallest neighbor of the broken cluster
    num = length(addedPins);
    LNBindex = [];
    NB = zeros(9,1);
    LNB = 2;
    for z = 1:num
        index = addedPins(z);
        i = ceil(index/max_rows);
        j = mod(index-1,max_col)+1;

        %Neighbors of pin                
        for m = 1:3
            for n = 1:3

                % Find the neighbor Value
                T = MLVandBL(i+(m-2),j+(n-2));

                %if neighbor value is a pin value, store as simple invader
                if T >= 2 || T <= 1
                    T = 2;
                end
                NB(n+(m-1)*3,1) = T;  
            end
        end

        %Check which neighbor is the smallest
        %and find the index                        
        for t = 1:9
            if LNB > NB(t,1)
                LNB = NB(t,1);
                if t < 4
                    m = mod(i-2, max_rows)+1;
                elseif t > 6
                    m = mod(i, max_rows)+1;
                else
                    m = i;
                end

                if mod(t,3) == 1
                    n = mod(j-2,max_col)+1;
                elseif mod(t,3) == 0
                    n = mod(j,max_col)+1;
                else
                    n = j;
                end

                LNBindex = (m-1)*max_rows+n;
            end
        end 
    end
end

