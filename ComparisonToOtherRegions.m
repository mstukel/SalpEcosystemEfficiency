clearvars
close all

load('CompOtherRegions.mat')
NPP = data.NPP_mgCM_2D_1_;
ProtistanGrazing = data.ProtistanGrazing_mgCM_2D_1_;
MesozooGrazing = data(:,6:10);
MesozooTP = data(:,13:17);
TunicateGrazing = data.GelatinousFilterFeederGrazing_mgCM_2D_1_;
TunicateSize = data.GelatinousFilterFeederMeanSize_mm_;
load('KeyParameters.mat')
SizeBinEdges = [ [0.5,1,2,5,10,20,50,100]/1000, ...
    [0.2,0.5,1,2,5,10,16,32,64,100,200,400,800,1600,3200]];


for cycle=1:19 %1:height(data)
    MicroGr = ProtistanGrazing(cycle);
    MesoGr = sum(table2array(MesozooGrazing(cycle,:)));
    SalpGr = TunicateGrazing(cycle);
    Zoo_TP = mean(table2array(MesozooTP(cycle,:)));
    Salp_TP = 2;  %Thaliaceans (Pyrosoma spinosum) were only present in Costa Rica Dome.  They were not substantially enriched in 15N relative to bulk seston.
    if isnan(data.ProtistanTrophicPosition(cycle))
        [Protist_TP,Protist2P] = CalcProtistTP(MicroGr,MesoGr,SalpGr,Zoo_TP,Salp_TP,ProtistGGE);
        Phy_TP = 1;
    else
        Protist_TP = data.ProtistanTrophicPosition(cycle);
        Phy_TP = data.PhototrophicFlagellateTrophicPosition(cycle);
    end
    
    Meso2P = zeros(1,length(SizeBinEdges)-1);
    MesoTP_track = ones(1,length(SizeBinEdges)-1)*2;
    for i=1:length(SizeBinEdges)-1
        if isnan(TunicateGrazing(cycle))==0
            if TunicateSize(cycle)>SizeBinEdges(i) & TunicateSize(cycle)<SizeBinEdges(i+1)
                SalpGr = TunicateGrazing(cycle);
                Salp_TP_track(1,i) = Salp_TP;
                Salp2P(1,i) = CalcZoo2PvariablePhy(SalpGr,Salp_TP,Protist_TP,Phy_TP,SalpGGE);
                if cycle < 6
                    Salp2P(1,i) = SalpGr*SalpGGE;  %Note that TP was very low for all pyrosomes in Costa Rica Dome
                end
            else
                Salp2P(1,i) = 0;
            end
        else
            Salp2P(1,i) = 0;
        end
        
        Mesozoomids = geomean([0.2,0.5,1,2,5;0.5,1,2,5,10]);
        for j=1:5
            if Mesozoomids(j)>SizeBinEdges(i) & Mesozoomids(j)<SizeBinEdges(i+1)
                MesoGr = table2array(MesozooGrazing(cycle,j));
                Meso_TP = table2array(MesozooTP(cycle,j));
                Meso_TP_track(1,i) = Meso_TP;
                Meso2P(1,i) = CalcZoo2PvariablePhy(MesoGr,Meso_TP,Protist_TP,Phy_TP,MesozooGGE);
                if Meso2P(1,i)<MesoGr*MesozooGGE
                    thatshouldnthappen
                end
            end
        end

    end

    [HTL2P, HTL_TP] = HigherTrophicLevels(Meso2P,Meso_TP_track,Salp2P,Salp_TP_track,SizeBinEdges,PredPreyRatioRange,HTLGGE);
    HTL2P_norm = sum(HTL2P)/NPP(cycle);

    HTL_track(cycle,:) = sum(HTL2P);

    load('SizeCutoff.mat')
    inds = find(SizeBinEdges(1:end-1)>=SizeCutoff);
    LargeProdNormOther(cycle) = sum(HTL2P_norm(inds));

end

metazoograzingfraction = (data.x0_2_0_5 + data.x0_5_1 + data.x1_2 + data.x2_5 + data.x_5 + data.GelatinousFilterFeederGrazing_mgCM_2D_1_)./data.NPP_mgCM_2D_1_;
LargeProdNormOther = [data.SurfaceChl_ugL_1_,TunicateGrazing./(TunicateGrazing + sum(table2array(MesozooGrazing)')'),LargeProdNormOther',data.NPP_mgCM_2D_1_,metazoograzingfraction,data.Picophytoplankton___,data.Picophytoplankton___]


load('LargeProdNorm.mat')
LargeProdNorm = [LargeProdNorm;LargeProdNormOther]
NZind = [1:5]; CRDind = [1:4]+5; EqPind = 5+5; NPSGind = 6+5; CCEind = [7:17]+5; GoMind = [18:19]+5;
fighandle = figure(82);
fighandle.Units = 'inches';
fighandle.Position = [5 3 5 3.8];
hold on
scatter(LargeProdNorm(NZind,1),LargeProdNorm(NZind,3),100,LargeProdNorm(NZind,2),'p','filled','MarkerEdgeColor','k')
scatter(LargeProdNorm(CRDind,1),LargeProdNorm(CRDind,3),80,LargeProdNorm(CRDind,2),'s','filled','MarkerEdgeColor','k')
scatter(LargeProdNorm(EqPind,1),LargeProdNorm(EqPind,3),50,LargeProdNorm(EqPind,2),'^','filled','MarkerEdgeColor','k')
scatter(LargeProdNorm(NPSGind,1),LargeProdNorm(NPSGind,3),50,LargeProdNorm(NPSGind,2),'v','filled','MarkerEdgeColor','k')
scatter(LargeProdNorm(CCEind,1),LargeProdNorm(CCEind,3),50,LargeProdNorm(CCEind,2),'o','filled','MarkerEdgeColor','k')
scatter(LargeProdNorm(GoMind,1),LargeProdNorm(GoMind,3),50,LargeProdNorm(GoMind,2),'d','filled','MarkerEdgeColor','k')
set(gca,'XScale','log')
h=colorbar
colormap(flipud(parula))
set(gca,'box','on')
legend('STF','CRD','EqP','NPSG','CCE','GoM')
xlabel('Surface Chlorophyll (\mug L^-^1)')
ylabel(['Ecosystem Transfer Efficiency',char(10),'(2^o Production of >',num2str(SizeCutoff),'-mm organisms / NPP)'])



fn = 'ComparisonToOtherRegions'
exportgraphics(gcf,[fn,'.pdf'],'Resolution',600)
exportgraphics(gcf,[fn,'.png'],'Resolution',600)


[rho,pval] = corr(LargeProdNorm(:,1),LargeProdNorm(:,3),'Type','Spearman')
fn = 'ManuscriptValues.txt'
fileID = fopen(fn,'a');
formatSpec = 'The Spearmans rank correlation between ecosystem transfer efficency and surface chlorophyll was %4.3f with p =  %4.3f.\n';
fprintf(fileID,formatSpec,[rho,pval])
fclose(fileID)


[rho,pval] = corr(LargeProdNorm(:,2),LargeProdNorm(:,3),'Type','Spearman')
fileID = fopen(fn,'a');
formatSpec = 'The correlation between ecosystem transfer efficency and the percentage of metazoan herbivory attributable to gelatinous filter feeders was %4.3f with p =  %4.3f.\n';
fprintf(fileID,formatSpec,[rho,pval])
fclose(fileID)

[rho,pval] = corr(LargeProdNorm(:,7),LargeProdNorm(:,3),'Type','Spearman')
fileID = fopen(fn,'a');
formatSpec = 'The correlation between ecosystem transfer efficency and the percentage of contribution of microphytoplankton to phyto biomass was %4.3f with p =  %4.3f.\n';
fprintf(fileID,formatSpec,[rho,pval])
fclose(fileID)

LargeProdNorm = array2table(LargeProdNorm,'VariableNames',{'SurfChl','GelPropHerbiv','EcosysTransferEfficiency','NPP','MetazooGraz_NPP','PicophyFrac','MicrophyFrac'});

Mdl = fitrgam(LargeProdNorm,'EcosysTransferEfficiency')