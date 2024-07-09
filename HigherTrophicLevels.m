function [HTL2P, HTL_TP] = HigherTrophicLevels(Meso2P,Meso_TP,Salp2P,Salp_TP,SizeBinEdges,PredPreyRatioRange,HTLGGE)

numTLs = 6;
SizeBinsMid = geomean([SizeBinEdges(1:end-1);SizeBinEdges(2:end)]);

HTL2P = zeros(numTLs,length(SizeBinsMid));
HTL_TP = zeros(numTLs,length(SizeBinsMid));

for i = 1:numTLs
    if i==1
        for j=1:length(SizeBinsMid)-1      %First walking over the size classes of metazoan (non-salp) zooplankton
            if Meso2P(j)>0
                HTL_TP(isnan(HTL_TP))=0;   %This is needed as a catch to fix divide by zeros for cells with no production
                PredRange = SizeBinsMid(j)*PredPreyRatioRange;   %Predator size range
                octavethickness0 = log2(PredPreyRatioRange(2))-log2(PredPreyRatioRange(1));  %This is the total number of octaves (factors of 2) in the predator range
                tmp = max(min(SizeBinEdges,PredRange(2)),PredRange(1));  %Dummy variable for setting up the proportion of secondary production that should be transferred
                octavethickness = log2(tmp(2:end))-log2(tmp(1:end-1));   %This is the octaves for size classes that correspond to potential predators of zooplankton j
                temp2P = Meso2P(j)*HTLGGE*octavethickness/octavethickness0;                        %This is the secondary production of the next trophic level, based on predation on zooplankton j
                temp_TP = Meso_TP(j)+1;                                   %This is the trophic position associated with predation on zooplankton j
                HTL_TP(i,:) = (HTL2P(i,:).*HTL_TP(i,:) + temp2P.*temp_TP) ./ (HTL2P(i,:) + temp2P);     %This combines with other zooplankton size classes to calculate the production-weighted trophic position 
                HTL2P(i,:) = HTL2P(i,:) + temp2P;
            end
        end
        for j=1:length(SizeBinsMid)-1      %Next walking over the size classes of salps
            if Salp2P(j)>0
                HTL_TP(isnan(HTL_TP))=0;   %This is needed as a catch to fix divide by zeros for cells with no production
                PredRange = SizeBinsMid(j)*PredPreyRatioRange;   %Predator size range
                octavethickness0 = log2(PredPreyRatioRange(2))-log2(PredPreyRatioRange(1));  %This is the total number of octaves (factors of 2) in the predator range
                tmp = max(min(SizeBinEdges,PredRange(2)),PredRange(1));  %Dummy variable for setting up the proportion of secondary production that should be transferred
                octavethickness = log2(tmp(2:end))-log2(tmp(1:end-1));   %This is the octaves for size classes that correspond to potential predators of zooplankton j
                temp2P = Salp2P(j)*HTLGGE*octavethickness/octavethickness0;                        %This is the secondary production of the next trophic level, based on predation on zooplankton j
                temp_TP = Salp_TP(j)+1;                                   %This is the trophic position associated with predation on zooplankton j
                HTL_TP(i,:) = (HTL2P(i,:).*HTL_TP(i,:) + temp2P.*temp_TP) ./ (HTL2P(i,:) + temp2P);     %This combines with other zooplankton size classes to calculate the production-weighted trophic position 
                HTL2P(i,:) = HTL2P(i,:) + temp2P;
            end
        end

    else  %Now walking over the size classes of higher trophic levels (i.e., HTLs eating other HTLs)
        for j=1:length(SizeBinsMid)-1
            if HTL2P(i-1,j)>0
                HTL_TP(isnan(HTL_TP))=0;   %This is needed as a catch to fix divide by zeros for cells with no production
                PredRange = SizeBinsMid(j)*PredPreyRatioRange;   %Predator size range
                octavethickness0 = log2(PredPreyRatioRange(2))-log2(PredPreyRatioRange(1));  %This is the total number of octaves (factors of 2) in the predator range
                tmp = max(min(SizeBinEdges,PredRange(2)),PredRange(1));  %Dummy variable for setting up the proportion of secondary production that should be transferred
                octavethickness = log2(tmp(2:end))-log2(tmp(1:end-1));   %This is the octaves for size classes that correspond to potential predators of zooplankton j
                temp2P = HTL2P(i-1,j)*HTLGGE*octavethickness/octavethickness0;                        %This is the secondary production of the next trophic level, based on predation on zooplankton j
                temp_TP = HTL_TP(i-1,j)+1;                                   %This is the trophic position associated with predation on zooplankton j
                HTL_TP(i,:) = (HTL2P(i,:).*HTL_TP(i,:) + temp2P.*temp_TP) ./ (HTL2P(i,:) + temp2P);     %This combines with other zooplankton size classes to calculate the production-weighted trophic position 
                HTL2P(i,:) = HTL2P(i,:) + temp2P;
            end
        end
    end
end

HTL_TP(isnan(HTL_TP))=0;   %This is needed as a catch to fix divide by zeros for cells with no production