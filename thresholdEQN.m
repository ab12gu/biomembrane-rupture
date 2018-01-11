function [ clusterSize ] = thresholdEQN( currentRadius, maxRadius )
%TRESE Calculates the treshold that determines the size of the cluster

clusterSize = 0.70*log2(0.5*(currentRadius/maxRadius)+1)+1;
        %2.3*log10(0.5*(CR/maxRadius)+1)+1;

end

