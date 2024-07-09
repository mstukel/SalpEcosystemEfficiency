clearvars
close all
warning('off','MATLAB:table:ModifiedAndSavedVarnames')

load('SecProd.mat')

SizeBinEdges = [ [0.5,1,2,5,10,20,50,100]/1000, ...
    [0.2,0.5,1,2,4,10,16,32,64,100,200,400,800,1600,3200]];


ProtistGGE = 0.3;
SalpGGE = 0.3;
MesozooGGE = 0.3;
HTLGGE = 0.3;
PredPreyRatioRange = [3,300];
origin = 'SecondaryProduction.m'
save('KeyParameters.mat','origin','ProtistGGE','SalpGGE','MesozooGGE','HTLGGE','PredPreyRatioRange')

load('SizeFracTP.mat')
load('SalpTP.mat')

FLOWCAM_HETEROTROPH_FRACTION = 0.33
Syn_C_cell = 255; %fgC / cell
Pro_C_cell = 36;  %fgC / cell


for cycle = 1:5
    cycle;
    NPP = GrowthGrazBalance.x14CPP(cycle);
    MicroGr = -GrowthGrazBalance.micrograz(cycle);
    MesoGr = -GrowthGrazBalance.zoop(cycle);
    SalpGr = -GrowthGrazBalance.salp(cycle);

    inds = find(SizeFracTPtable.Cycle==cycle);
    Zoo_TP = mean(SizeFracTPtable.TP(inds));

    inds = find(SalpTPtable.Cycle==cycle);
    Salp_TP = nanmean(SalpTPtable.TP(inds));
    if Salp_TP == 0 | isnan(Salp_TP)
        Salp_TP = 2.3;  %This doesnt matter, because SalpGr was 0 anytime Salp_TP was NaN;
    end

    [Protist_TP,Protist2P] = CalcProtistTP(MicroGr,MesoGr,SalpGr,Zoo_TP,Salp_TP,ProtistGGE);

    %Salp Grazing Rate Calculation
    SalpSizeBinMid = Field.SizeBinMid;
    if cycle == 1
        SalpCommunityGutPig = Field.C1_CommunityGutPig;
        ProportionalSalpGutPig = SalpCommunityGutPig/sum(SalpCommunityGutPig);
    elseif cycle == 2
        SalpCommunityGutPig = Field.C2_CommunityGutPig;
        ProportionalSalpGutPig = SalpCommunityGutPig/sum(SalpCommunityGutPig);
    elseif cycle == 4
        SalpCommunityGutPig = Field.C4_CommunityGutPig;
        ProportionalSalpGutPig = SalpCommunityGutPig/sum(SalpCommunityGutPig);
    else
        ProportionalSalpGutPig = zeros(size(SalpSizeBinEdges));  %Cycles with no salps
    end
    for i=1:length(SizeBinEdges)-1
        inds = find(SalpSizeBinMid>SizeBinEdges(i) & SalpSizeBinMid<SizeBinEdges(i+1));
        SalpGr_Size(1,i) = sum(ProportionalSalpGutPig(inds))*SalpGr;
    end

    %Salp Secondary Production Calculations
    if cycle==1 | cycle==2 | cycle==4
        for i=1:length(SizeBinEdges)-1
            inds = find(SalpTPtable.Cycle==cycle & SalpTPtable.Length>SizeBinEdges(i) & SalpTPtable.Length<SizeBinEdges(i+1));
            Salp_TP(i,1) = nanmean(SalpTPtable.TP(inds));
            if Salp_TP(i,1) == 0 | isnan(Salp_TP(i,1))
                Salp_TP(i,1) = nanmean(SalpTPtable.TP(SalpTPtable.Cycle==cycle));  %If no salps in that size class for that cycle, replacing with the mean over that cycle
            end

            Salp2P(1,i) = CalcZoo2P(SalpGr_Size(1,i),Salp_TP(i,1),Protist_TP,SalpGGE);

        end

        SalpTotalIngestion = sum(Salp2P)/SalpGGE;
        SalpFractionProtist(cycle)=1 - sum(SalpGr_Size)/SalpTotalIngestion;   %Calculating the fraction of protistan zooplankton in salp diets
    else
        Salp2P = zeros(1,length(SizeBinEdges)-1);  %Non-salp cycles
        SalpFractionProtist(cycle) = NaN;
    end


    %Size-fractionated zooplankton grazing rates
    ind = find(zoopgrazdata.Cycle==cycle);
    MEAN = mean([zoopgrazdata.x0_2_0_5mm_4(ind), zoopgrazdata.x0_5_1mm_4(ind), zoopgrazdata.x1_2mm_4(ind), zoopgrazdata.x2_4mm_4(ind), zoopgrazdata.x_4mm_4(ind)]);
    for i=1:length(MEAN)
        ProportionalZoopGutPig(i) = MEAN(i)/sum(MEAN);
    end
    SizeFracSizeBinMid = geomean([0.2 0.5 1 2 5; 0.5 1 2 5 10]);
    for i=1:length(SizeBinEdges)-1
        inds = find( SizeFracSizeBinMid>SizeBinEdges(i) &  SizeFracSizeBinMid<SizeBinEdges(i+1));
        MesoGr_Size(1,i) = sum(ProportionalZoopGutPig(inds))*MesoGr;
    end


    %Size-fractionated zooplankton secondary Production Calculations
    for i=1:length(SizeBinEdges)-1
        inds = find(SizeFracTPtable.Cycle==cycle & SizeFracTPtable.MedSize>SizeBinEdges(i) & SizeFracTPtable.MedSize<SizeBinEdges(i+1));
        Meso_TP(i,1) = nanmean(SizeFracTPtable.TP(inds));
        if Meso_TP(i,1) == 0 | isnan(Meso_TP(i,1))
            Meso_TP(i,1) = 2.5;  %Just setting a typical number which applies to size classes for which we did not sample, although such size classes will have no grazing, so it won't matter (I just need a non-zero / non-NaN number so that it doesn't blow up the calculations)
        end
        Meso2P(1,i) =CalcZoo2P(MesoGr_Size(1,i),Meso_TP(i,1),Protist_TP,MesozooGGE);
    end

    MesoTotalIngestion = sum(Meso2P)/MesozooGGE;
    MesoFractionProtist(cycle)=1 - sum(MesoGr_Size)/MesoTotalIngestion;   %Calculating the fraction of protistan zooplankton in size-fractionated zooplankton diets

    %Secondary Production of Higher Trophic Levels Calculations
    [HTL2P, HTL_TP] = HigherTrophicLevels(Meso2P,Meso_TP,Salp2P,Salp_TP,SizeBinEdges,PredPreyRatioRange,HTLGGE);
    clear i ind inds MEAN MesoGr MicroGr MesoGr_Size ProportionalSalpGutPig ProportionalZoopGutPig SalpCommunityGutPig SalpGr SalpGr_Size SizeFracSizeBinMid tmp Zoo_TP

    %Picoeukaryotic phytoplankton size and abundance
    ind = find(NBSStable_keep.Cycle)==cycle;
    NBSStable = NBSStable_keep(ind,:);
    picoeuk_biomass = zeros(1,length(SizeBinEdges)-1);
    for i=1:length(SizeBinEdges)-1
        if SizeBinEdges(i)<4*10^-3  %For the 4 - 8 um size ranges on up, flow cam is more reliable
            for j=6:width(NBSStable)
                var = NBSStable.Properties.VariableNames(j);
                midpt = str2num(var{1});
                if midpt>SizeBinEdges(i)*1000 & midpt<SizeBinEdges(i+1)*1000  %Times 1000 is converting from microns to millimeters
                    nbss = mean(table2array(NBSStable(:,j)));
                    binwid = midpt*scaling - midpt/scaling;
                    picoeuk_biomass(i) = nbss*binwid;  %picograms C / mL
                end
            end
        end
    end

    %Cyanobacteria Biomass
    cyano_biomass = zeros(1,length(SizeBinEdges)-1);
    ind = find(flowcyt_cyano.CYCLE==cycle);
    pro_abund = flowcyt_cyano.PRO_ML(ind);
    pro_abund(find(isnan(pro_abund)))=0;  %Missing values were Prochlorococcus not detected
    syn_abund = flowcyt_cyano.SYN_ML(ind);
    pro_abund = mean(pro_abund);  %cells/mL
    syn_abund = mean(syn_abund);  %cells/mL
    pro_biom = pro_abund*Pro_C_cell/1000;  %Converting femto to pico: final units are picograms C / mL
    syn_biom = syn_abund*Syn_C_cell/1000;  %Converting femto to pico: final units are picograms C / mL
    if SizeBinEdges(1)==0.5/1000
        cyano_biomass(1) = pro_biom;  %Putting Prochlo in the 0.5 - 1 um size bin
        cyano_biomass(2) = syn_biom;  %Putting Syn in the 1 - 2 um size bin
    end

    %FlowCam Biomasses;
    ind = find(strcmp(flowcam.Cycle,['Cycle ',num2str(cycle)]));
    flowcam_other = zeros(1,length(SizeBinEdges)-1);
    flowcam_centric = zeros(1,length(SizeBinEdges)-1);
    flowcam_pennate = zeros(1,length(SizeBinEdges)-1);
    flowcam_dino = zeros(1,length(SizeBinEdges)-1);
    flowcam_silico = zeros(1,length(SizeBinEdges)-1);
    flowcam_ciliate = zeros(1,length(SizeBinEdges)-1);
    for i=1:length(SizeBinEdges)-1
        for j=56:width(flowcam)  %Skippping to the abundance columns
            var_name = flowcam.Properties.VariableNames(j);
            var_name = var_name{1};
            ISNUMBER = [];
            for k=1:length(var_name)
                if str2num(var_name(k))>-1 & isreal(str2num(var_name(k)))
                    ISNUMBER(k)=1;
                else
                    ISNUMBER(k)=0;
                end
            end
            tmp = find(ISNUMBER==1);
            start1=tmp(1);
            end1=tmp(end);
%             for k=2:length(tmp)
%                 if tmp(k)-1~=tmp(k-1)
%                     end1=tmp(k-1);
%                     %start2=tmp(k);
%                     %end2=tmp(end);
%                     break
%                 end
%             end
            Group = var_name(1:start1-1);
            lower_lim = str2num(var_name(start1:end1));
            %upper_lim = str2num(var_name(start2:end2));
            upper_lim = lower_lim*2;
            midpt = geomean([lower_lim,upper_lim])/1000;
            binwid = upper_lim - lower_lim;
            if SizeBinEdges(i)<midpt & SizeBinEdges(i+1)>=midpt & midpt>4/1000 & strcmp(Group,'Other_')    %Note that the 4 is in there, because we only trust flowcam for >4-um cells, below that we rely on the picoeuk estimates from flow cytometry
                flowcam_other(i) = flowcam_other(i)+table2array(flowcam(ind,j));
            elseif SizeBinEdges(i)<midpt & SizeBinEdges(i+1)>=midpt & midpt>4/1000 & strcmp(Group,'Pennate_')    %Note that the 4 is in there, because we only trust flowcam for >4-um cells, below that we rely on the picoeuk estimates from flow cytometry
                flowcam_pennate(i) = flowcam_other(i)+table2array(flowcam(ind,j));
            elseif SizeBinEdges(i)<midpt & SizeBinEdges(i+1)>=midpt & midpt>4/1000 & strcmp(Group,'Centric_')    %Note that the 4 is in there, because we only trust flowcam for >4-um cells, below that we rely on the picoeuk estimates from flow cytometry
                flowcam_centric(i) = flowcam_other(i)+table2array(flowcam(ind,j));
            elseif SizeBinEdges(i)<midpt & SizeBinEdges(i+1)>=midpt & midpt>4/1000 & strcmp(Group,'Silico_')    %Note that the 4 is in there, because we only trust flowcam for >4-um cells, below that we rely on the picoeuk estimates from flow cytometry
                flowcam_silico(i) = flowcam_other(i)+table2array(flowcam(ind,j));
            elseif SizeBinEdges(i)<midpt & SizeBinEdges(i+1)>=midpt & midpt>4/1000 & strcmp(Group,'Dino_')    %Note that the 4 is in there, because we only trust flowcam for >4-um cells, below that we rely on the picoeuk estimates from flow cytometry
                flowcam_dino(i) = flowcam_other(i)+table2array(flowcam(ind,j));
            elseif SizeBinEdges(i)<midpt & SizeBinEdges(i+1)>=midpt & midpt>4/1000 & strcmp(Group,'Ciliate_')    %Note that the 4 is in there, because we only trust flowcam for >4-um cells, below that we rely on the picoeuk estimates from flow cytometry
                flowcam_ciliate(i) = flowcam_other(i)+table2array(flowcam(ind,j));
            end

        end
    end

        flowcam_phybiomass = flowcam_centric + flowcam_pennate + flowcam_silico + flowcam_dino*(1-FLOWCAM_HETEROTROPH_FRACTION) + flowcam_other*(1-FLOWCAM_HETEROTROPH_FRACTION);
        flowcam_hetbiomass = flowcam_ciliate + flowcam_dino*FLOWCAM_HETEROTROPH_FRACTION + flowcam_other*FLOWCAM_HETEROTROPH_FRACTION;
    phy_biomass = cyano_biomass + picoeuk_biomass + flowcam_phybiomass;
    Pico_fraction(cycle,1) = sum(phy_biomass(1:2))./sum(phy_biomass);
    Micro_fraction(cycle,1) = sum(phy_biomass(6:end))./sum(phy_biomass);

    PHYTOPLANKTON_BIOMASS(cycle,1) = sum(phy_biomass);

    hetpro_biomass = flowcam_hetbiomass;
    phy2P = NPP*phy_biomass/sum(phy_biomass);
    hetpro2P = Protist2P*hetpro_biomass/sum(hetpro_biomass);
    HETEROTROPHICPROTIST_BIOMASS(cycle,1) = sum(hetpro_biomass);

    Secondary_Production(cycle,:) = hetpro2P + Meso2P + Salp2P + sum(HTL2P);
    Total_Production(cycle,:) = phy2P + hetpro2P + Meso2P + Salp2P + sum(HTL2P);
    Secondary_Production_norm(cycle,:) = Secondary_Production(cycle,:)/NPP;
    Total_Production_norm(cycle,:) = Total_Production(cycle,:)/NPP;

    Meso2P_track(cycle,1) = sum(Meso2P);
    Protist2P_track(cycle,1) = sum(Protist2P);
    Salp2P_track(cycle,1) = sum(Salp2P);
end

origin = 'SecondaryProduction.m'
save('PicoMicroPhyFractions.mat','Pico_fraction','Micro_fraction','origin')
units = 'pgC/mL';
save('ProtistBiomass.mat','PHYTOPLANKTON_BIOMASS','HETEROTROPHICPROTIST_BIOMASS','origin','units')
clear units PHYTOPLANKTON_BIOMASS HETEROTROPHICPROTIST_BIOMASS
fn = 'ManuscriptValues.txt'

fileID = fopen(fn,'a');
formatSpec = 'The proportion of protistan zooplankton in salp guts was %4.4f, %4.4f, and %4.4f for Cycles 1, 2, and 4, respectively.\n';
fprintf(fileID,formatSpec,[SalpFractionProtist(1),SalpFractionProtist(2),SalpFractionProtist(4)])
fclose(fileID)

fileID = fopen(fn,'a');
formatSpec = 'The proportion of protistan zooplankton in size-fractionated zooplankton guts ranged from %4.4f to %4.4f .\n';
fprintf(fileID,formatSpec,[min(MesoFractionProtist),max(MesoFractionProtist)])
fclose(fileID)

fileID = fopen(fn,'a');
formatSpec = 'The proportion of protistan zooplankton in salp guts averaged %4.4f +/- %4.4f  when salps were present.\n';
fprintf(fileID,formatSpec,[nanmean(SalpFractionProtist),nanstd(SalpFractionProtist)])
fclose(fileID)

fileID = fopen(fn,'a');
formatSpec = 'The proportion of protistan zooplankton in size-fractionated zooplankton averaged %4.4f +/- %4.4f .\n';
fprintf(fileID,formatSpec,[mean(MesoFractionProtist),std(MesoFractionProtist)])
fclose(fileID)


fileID = fopen(fn,'a');
formatSpec = 'Protistan secondary production ranged from %4.4f to %4.4f mg C m^-^2 d^-^1 .\n';
fprintf(fileID,formatSpec,[min(Protist2P_track),max(Protist2P_track)])
fclose(fileID)

fileID = fopen(fn,'a');
formatSpec = 'Mesozooplankton secondary production ranged from %4.4f to %4.4f mg C m^-^2 d^-^1 .\n';
fprintf(fileID,formatSpec,[min(Meso2P_track),max(Meso2P_track)])
fclose(fileID)

fileID = fopen(fn,'a');
formatSpec = 'Salp secondary production ranged from %4.4f to %4.4f mg C m^-^2 d^-^1 .\n';
fprintf(fileID,formatSpec,[min(Salp2P_track([1,2,4])),max(Salp2P_track([1,2,4]))])
fclose(fileID)



clearvars -except SizeBinEdges Secondary_Production Total_Production Secondary_Production_norm Total_Production_norm GrowthGrazBalance
SizeBinMid = geomean([SizeBinEdges(1:end-1);SizeBinEdges(2:end)]);
Octaves = (SizeBinEdges(2:end)./SizeBinEdges(1:end-1))/2;  %This is the number of octaves (i.e., factors of 2)
origin = 'SecondaryProduction.m'
save('SecondaryProduction.mat','origin','SizeBinMid','Octaves','Secondary_Production_norm','Secondary_Production','Total_Production')

%---------------------------------------------------------------------------------------------------------------------
%-------------Secondary Production Figure-----------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------


load('CycleColors.mat')

fighandle = figure(172);
fighandle.Units = 'inches';
fighandle.Position = [3 1 3.8 5];

subplot(2,1,1)
hold on
for i=1:5
    if i==1 | i==2 | i==4
        plot(SizeBinMid,Total_Production(i,:)./Octaves,'-ok','Color',cols(i,:),'MarkerFaceColor',cols(i,:),'MarkerEdgeColor','k')
    else
        plot(SizeBinMid,Total_Production(i,:)./Octaves,'-dk','Color',cols(i,:),'MarkerFaceColor',cols(i,:),'MarkerEdgeColor',[0.5 0.5 0.5])
    end
end
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'box','on')
% h=legend('C1 (salp)','C2 (salp)','C3','C4 (salp)','C5','Location','SouthWest');
% set(h,'FontSize',8)
set(gca,'FontSize',9)
xlabel('Organism Size (mm)')
ylabel(['Biomass Production',char(10), '(mg C m^-^2 d^-^1 / octave)'])
xlim([min(SizeBinEdges),max(SizeBinEdges)])
text(SizeBinEdges(end)*0.4,max(max(Total_Production)./Octaves)*0.6,'a')
set(gca,'YTick',10.^[-1 0 1 2])
set(gca,'XTick',10.^[-4 -3 -2 -1 0 1 2 3 4])

subplot(2,1,2)
hold on
for i=1:5
    if i==1 | i==2 | i==4
        plot(SizeBinMid,Secondary_Production_norm(i,:)./Octaves,'-ok','Color',cols(i,:),'MarkerFaceColor',cols(i,:),'MarkerEdgeColor','k')
    else
        plot(SizeBinMid,Secondary_Production_norm(i,:)./Octaves,'-dk','Color',cols(i,:),'MarkerFaceColor',cols(i,:),'MarkerEdgeColor',[0.5 0.5 0.5])
    end
end
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'box','on')
h=legend('C1 (salp)','C2 (salp)','C3','C4 (salp)','C5','Location','SouthWest');
set(h,'FontSize',7)
set(gca,'FontSize',9)
ylabel(['Secondary Production',char(10),'Normalized to NPP',char(10),'(per octave)'])
xlim([min(SizeBinEdges),max(SizeBinEdges)])
xlabel('Organism Size (mm)')
text(SizeBinEdges(end)*0.4,max(max(Secondary_Production_norm)./Octaves)*0.7,'b')
set(gca,'YTick',10.^[-4 -3 -2 -1])
set(gca,'XTick',10.^[-4 -3 -2 -1 0 1 2 3 4])
set(gca,'YLim',[10^-4 0.4])

fn = 'SecondaryProduction.TotProd&NormalizedSecProd'
exportgraphics(gcf,[fn,'.pdf'],'Resolution',600)
exportgraphics(gcf,[fn,'.png'],'Resolution',600)



%---------------------------------------------------------------------------------------------------------------------
%-------------Manuscript Text -----------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------

SalpCycles = [1,2,4];
NonSalpCycles = [3,5];
SizeCutoff = 200;  %mm
save('SizeCutoff.mat','SizeCutoff','origin')

inds = find(SizeBinEdges(1:end-1)>=SizeCutoff);
LargeProdNorm = sum(Total_Production_norm(:,inds)')';

fn = 'ManuscriptValues.txt'
fileID = fopen(fn,'a');
formatSpec = 'The ratio of secondary production by > %4.0f mm individuals to total primary production was %4.3f , %4.3f , and %4.3f for salp cycles.\n';
fprintf(fileID,formatSpec,[SizeCutoff,LargeProdNorm(SalpCycles(1)),LargeProdNorm(SalpCycles(2)),LargeProdNorm(SalpCycles(3))])
fclose(fileID)

fileID = fopen(fn,'a');
formatSpec = 'For non-salp cycles it was was %4.4f  and %4.4f .\n';
fprintf(fileID,formatSpec,[LargeProdNorm(NonSalpCycles(1)),LargeProdNorm(NonSalpCycles(2))])
fclose(fileID)


    MesoGr = -GrowthGrazBalance.zoop;
    SalpGr = -GrowthGrazBalance.salp;
    PerGrSalp = SalpGr./(SalpGr+MesoGr);
SurfChl = [0.88; 0.41; 2.21; 1.3; 0.21];  %Fluorometric chlorophyll from Yingling et al. (in prep).
load('PicoMicroPhyFractions.mat')
LargeProdNorm = [SurfChl, PerGrSalp, LargeProdNorm, GrowthGrazBalance.x14CPP, (-GrowthGrazBalance.zoop-GrowthGrazBalance.salp)./GrowthGrazBalance.x14CPP, Pico_fraction, Micro_fraction];
save('LargeProdNorm.mat','LargeProdNorm','origin')















