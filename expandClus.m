function [pins, circle] = expandClus(pinB, pins, circle, LNBindex, thres, max_rows, max_col)
%Expand the cluster from the broken point
    Len = length(pins)+1;
    FF = find(LNBindex == pins(1,:));

    pins(1,Len) = LNBindex;

    i = ceil(LNBindex/max_rows);
    j = mod(LNBindex-1,max_col)+1;
    circle(i, j) = 1;

    %frame=getframe(gcf);
    
    for m = 1:3
        for n = 1:3
            NBi = i + (m-2);
            NBj = j + (n-2);

            T = circle(NBi, NBj);

            if T < thres && T > 1
                index = (NBi-1)*max_rows+NBj;
                LNBindex = index;
                [pins, circle] = expandClus(pinB, pins, circle, LNBindex, thres, max_rows, max_col);
            end
        end
    end
end

