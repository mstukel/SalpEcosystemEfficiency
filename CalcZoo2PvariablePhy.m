function Zoo2P = CalcZoo2PvariablePhy(GutPig,Zoo_TP,Protist_TP,Phy_TP,GGE)

DF_Zoo_Phy = (Zoo_TP - Protist_TP - 1)/(Phy_TP - Protist_TP);
%DF_Zoo_HetPro = 1 - DF_Zoo_Phy;

Zoo2P = GutPig/DF_Zoo_Phy*GGE;