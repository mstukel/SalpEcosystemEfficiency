function [Protist_TP,Protist2P,ProductionWeightedProtistTP ] = CalcProtistTP(MicroGr,MesoGr,SalpGr,Zoo_TP,Salp_TP,ProtistGGE)

%-------------------------------------------------------------------------
% This function calculates protist trophic position (Protist_TP) and
% protist secondary production (Protist2P) from the carbon-based grazing
% rates of microzooplankton (MicroGr), mesozoopankton (MesoGr), and salps
% (SalpGr); the gross growth efficiency of protists (ProtistGGE), and the
% trophic positions of mesozooplankton (Zoo_TP) and salps (Salp_TP).
% Equations are explained in the "Calculating production as a function of 
% size and trophic position" section of Stukel et al.
%-------------------------------------------------------------------------


% clearvars
% close all
% MicroGr = 370.8363;
% MesoGr = 114.1297;
% SalpGr = 267.4123;
% Salp_TP = 2.4903;
% Zoo_TP = 2.5171;
% ProtistGGE = 0.3;


Protist_TP = 100;  %Initial Guess
Phy_TP = 1;
%Protist2P = zeros(MaxNumProtistTrophicSteps,1);
Protist2P(1) = MicroGr*ProtistGGE;  %Initial Guess is that all production is at the second trophic level


for iter = 1:10

    %Calculating dietary fractions and total ingestion for non-salp metazoan zooplankton
    DF_Zoo_Phy = (Zoo_TP - Protist_TP - 1)/(Phy_TP - Protist_TP);
    DF_Zoo_HetPro = 1 - DF_Zoo_Phy;
    DF_Zoo_Phy = min([DF_Zoo_Phy,1]);
    DF_Zoo_Phy = max([DF_Zoo_Phy,0]);
    Zoo_TotalIngestion = MesoGr./DF_Zoo_Phy;

    %Calculating dietary fractions and total ingestion for non-salp metazoan zooplankton
    DF_Salp_Phy = (Salp_TP - Protist_TP - 1)/(Phy_TP - Protist_TP);
    DF_Salp_HetPro = 1 - DF_Salp_Phy;
    DF_Salp_Phy = min([DF_Salp_Phy,1]);
    DF_Salp_Phy = max([DF_Salp_Phy,0]);
    Salp_TotalIngestion = SalpGr./DF_Salp_Phy;

    %Calculating the amount of energy that must be dissipated through protistan respiration, excretion, and defecation
    Dissipation_Protist = MicroGr - Salp_TotalIngestion*DF_Salp_HetPro - Zoo_TotalIngestion*DF_Zoo_HetPro;

    %Calculating protistan trophic position, given the amount of energy they must dissipate and the phytoplankton carbon they consumed (and their GGE)
    Protist_TP = log( -(Dissipation_Protist-MicroGr)/MicroGr ) / log(ProtistGGE) + 1;   %This is a solution to the equation: Dissipation = (1-GGE)*MicroGr*(1-GGE^(TP-1))/(1-GGE)
    
    %This is a catch for issues that occur when total carbon consumption rates cannot be satisfied by the measured protistan grazing, when that occurs, we assume that the secondary production of protistan herbivores must be completely consumed by metazoan zooplankton
    if Protist_TP<2
        Protist_TP=2;
    end

    %For the secondary production v. trophic level plots, I will want to know the production-weighted average trophic level of the heterotrophic protists.  Calculating that here.
    Protist2P = ProtistGGE * MicroGr * ( 1- ProtistGGE^(Protist_TP-1))/(1-ProtistGGE);
    for i = 1:(Protist_TP - 2)*100+1
        tp(i) = (i-1)/100+2;
        if i==1
            Protist2P_temp(i) = ProtistGGE * MicroGr * ( 1- ProtistGGE^(tp(i)-1))/(1-ProtistGGE);
        else
            Protist2P_temp(i) = ProtistGGE * MicroGr * ( 1- ProtistGGE^(tp(i)-1))/(1-ProtistGGE) - sum(Protist2P_temp(1:i-1));
        end
    end
    ProductionWeightedProtistTP = sum(Protist2P_temp.*tp)/sum(Protist2P_temp);

    %

end
            



