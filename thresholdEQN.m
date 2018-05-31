function [ clusterSize ] = thresholdEQN( pin_distance_to_center, maxRadius )
%TRESE Calculates the treshold that determines the size of the cluster

clusterSize = 0.70*log2(0.5*(pin_distance_to_center/maxRadius)+1)+1;
        %2.3*log10(0.5*(CR/maxRadius)+1)+1;

end

