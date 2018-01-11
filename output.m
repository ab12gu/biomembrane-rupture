function [] = output(wholeCirc, pinProb, clusProb, maxRadius, frameSize, ...
    expRadius, videoSave, pins, MLVnodes)
%OUTPUT Saves a video of the output array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Show All Pins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pinLength = size(pins,2);
    
    for node = 1:pinLength
        index = pins(1,node);
        i = ceil(index/frameSize);
        j = mod(index-1,frameSize)+1;
        if pins(8,node) > 0
            if pins(10,node) ~= 1
                wholeCirc(i,j) = 3; %If pin (clustering)
            else
                wholeCirc(i,j) = 4; %If fractured node
            end
        else
            wholeCirc(i,j) = 3; %If only pin site
        end
    end

    for node = 1:length(MLVnodes)
        index = MLVnodes(1,node);
        i = ceil(index/frameSize);
        j = mod(index-1,frameSize)+1;
        wholeCirc(i,j) = 1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visual Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    image(ceil(wholeCirc)+1)
    axis image

    title({['Pin Density: ',  num2str(pinProb*100), '%'];...
        ['Cluster Density: ', num2str(clusProb*100), '%' ]});
    set(gca,'xtick',[]);
    set(gca,'YTick',[]);
    xlabel({['Maxium Radius:' num2str(maxRadius), ' || Frame Size:' ...
        num2str(frameSize), ':', num2str(frameSize)];...
        [' Radius of Displaced Bi-layer: ', ...
        num2str(round(expRadius))]});

    %Frame Settings
    width= 850;
    height= 600;
    set(gcf,'position',[300,50,width,height])

    %Colorbar Settings:
    h = colorbar;  
     set(h, 'XTick', [1.5, 2.5, 3.3, 3.7, 4.5, 5.5])
     set(h,'XTickLabel', {   'Background',...
                             'Source (MLV)',...
                             'Bi-Layer',...
                             'Displaced',...
                             'Pinning Sites',...
                             'Fractured Area'})
    title(h, '                    Colormap Legend')

    frame=getframe(gcf);
    
    if videoSave == 1
        writeVideo(v,frame);
    end
end

