addpath '..\..\..\..\Misc Oceanography\Matlab Functions'

clearvars
close all
warning('off','MATLAB:table:ModifiedAndSavedVarnames')
load('TP.mat')

for i=1:length(SrcAA)
    SizeFrac_Src(:,i) = eval([char('SizeFrac.'),char(SrcAA(i))]);
end
SizeFrac_Src = nanmean(SizeFrac_Src')';

for i=1:length(TrAA)
    SizeFrac_Tr(:,i) = eval([char('SizeFrac.'),char(TrAA(i))]);
end
SizeFrac_Tr = nanmean(SizeFrac_Tr')';
SizeFrac_TP = (SizeFrac_Tr - SizeFrac_Src - Beta)/TDF_eco + 1;


SizeFracTPtable = [SizeFrac.Cycle,SizeFrac.MedSize,SizeFrac_TP];
SizeFracTPtable(find(isnan(SizeFracTPtable(:,3))),:)=[];
SizeFracTPtable = array2table(SizeFracTPtable,'VariableNames',{'Cycle','MedSize','TP'});
origin = 'TrophicPositions.m';
save('SizeFracTP.mat','SizeFracTPtable','origin')



for i=1:length(SrcAA)
    Salps_Src(:,i) = str2double(eval([char('Salps.'),char(SrcAA(i))]));
end
Salps_Src = nanmean(double(Salps_Src)')';

for i=1:length(TrAA)
    Salps_Tr(:,i) = str2double(eval([char('Salps.'),char(TrAA(i))]));
end
Salps_Tr = nanmean(double(Salps_Tr)')';
Salps_TP = (Salps_Tr - Salps_Src - Beta - TDF_salp)/TDF_eco + 2;

SalpTPtable = [Salps.Cycle,Salps.Length,Salps_TP];
SalpTPtable = array2table(SalpTPtable,'VariableNames',{'Cycle','Length','TP'});
SalpTPtable = [SalpTPtable,Salps(:,2:4)];
inds = find(SalpTPtable.Body1Gut0Hyp2==1);
SalpTPtable = SalpTPtable(inds,:);
save('SalpTP.mat','SalpTPtable','origin')

%---------------------------------------------------------------------
%-------------Protistan Trophic Steps---------------------------------
%---------------------------------------------------------------------


TP_ala = (SizeFrac.Ala - SizeFrac.Phe - Beta_ala)/TEF_ala + 1;
TP_glx = (SizeFrac.Glx - SizeFrac.Phe - Beta_glx)/TEF_glx + 1;

SizeFrac_deltaTP = TP_ala - TP_glx;


SizeFrac_deltaTP(find(isnan(SizeFrac_deltaTP)))=[];
SizeFrac_deltaTP_mean = mean(SizeFrac_deltaTP);
SizeFrac_deltaTP_se = std(SizeFrac_deltaTP)/sqrt(length(SizeFrac_deltaTP));

fn = 'ManuscriptValues.txt'
fileID = fopen(fn,'a');
%fprintf(fileID,['Calculated Trophic Discrimination Factor for salps was ',num2str(TDF_salp,2),' +/- ',num2str(TDF_salp_se,2)])
formatSpec = 'The mean number of protistan trophic steps calculated from size-fractionated mesozooplankton samples averaged  %4.2f +/- %4.2f \n';
fprintf(fileID,formatSpec,[SizeFrac_deltaTP_mean,SizeFrac_deltaTP_se])
fclose(fileID)


%---------------------------------------------------------------------
%-------------SizeFractionated and Salps Figure---------------------------------
%---------------------------------------------------------------------


fighandle = figure(1132);
fighandle.Units = 'inches';
fighandle.Position = [3 1 8 4];

subplot(1,2,1)
hold on
REORG = NaN(50,5);
for i=1:5
    inds = find(SizeFrac.Cycle==i);
    temp = SizeFrac_TP(inds);
    REORG(1:length(temp),i)=temp;
end
boxplot(REORG)
set(gca,'box','on')
ylabel('Zooplankton Trophic Position (TP_A_A)')
xlabel('Lagrangian Experiment')
set(gca,'XTickLabel',{'C1 - SA-Sc','C2 - SA','C3 - ST','C4 - ST','C5 - SA'})
ylim([1.3 3.8])
set(gca,'FontSize',9)
text(5.2,1.4,'a','FontSize',10)


subplot(1,2,2)
hold on
for i=1:6
    if i==6
        inds = find(Salps.Cycle==1.5 & Salps.Body1Gut0Hyp2==1 & strcmp(Salps.Species,'ST'));  %Salps collected outside of cycles
        plot(Salps.Length(inds),Salps_TP(inds),'ok','MarkerFaceColor',[0.7 0.7 0.7])
        inds = find(Salps.Cycle==1.5 & Salps.Body1Gut0Hyp2==1 & strcmp(Salps.Species,'TV'));  %Salps collected outside of cycles
        plot(Salps.Length(inds),Salps_TP(inds),'dk','MarkerFaceColor',[0.7 0.7 0.7])
        inds = find(Salps.Cycle==1.5 & Salps.Body1Gut0Hyp2==1 & strcmp(Salps.Species,'PC'));  %Salps collected outside of cycles
        plot(Salps.Length(inds),Salps_TP(inds),'sk','MarkerFaceColor',[0.7 0.7 0.7])
    else
        inds = find(Salps.Cycle==i & Salps.Body1Gut0Hyp2==1 & strcmp(Salps.Species,'ST'));
        plot(Salps.Length(inds),Salps_TP(inds),'ok','MarkerFaceColor',cols(i,:))
        inds = find(Salps.Cycle==i & Salps.Body1Gut0Hyp2==1 & strcmp(Salps.Species,'PC'));
        plot(Salps.Length(inds),Salps_TP(inds),'sk','MarkerFaceColor',cols(i,:))
        inds = find(Salps.Cycle==i & Salps.Body1Gut0Hyp2==1 & strcmp(Salps.Species,'SZ'));
        plot(Salps.Length(inds),Salps_TP(inds),'^k','MarkerFaceColor',cols(i,:))
    end
end
xlabel('Length (mm)')
ylabel('Salp Trophic Position (TP_A_A)')
set(gca,'box','on')
Bounds = [50 147; 2.9 3.76];
Shapes = {'o','o','o','o','o','o'};
Labels = {'C1';'C2';'C3';'C4';'C5';'other'};
[output] = MakeLegend(Bounds,Shapes,flipud([cols;0.7,0.7,0.7]),flipud(Labels),8);
Bounds(1,1) = 90
Shapes = {'o','d','s','^'};
Labels = {'\itS. thompsoni';'\itT. vagina';'\itP. confoederata';'\itS. zonaria'};
[output] = MakeLegend(Bounds,fliplr(Shapes),flipud(cols(1:4,:)*0+1),flipud(Labels),8);
text(143,1.4,'b','FontSize',10)
set(gca,'FontSize',9)
ylim([1.3 3.8])


fn = 'TrophicPositions.SizeFrac&Salps'
exportgraphics(gcf,[fn,'.pdf'],'Resolution',600)
exportgraphics(gcf,[fn,'.png'],'Resolution',600)









%------------------------------------------------------------------------
%----------Values for Manuscript-----------------------------------------
%------------------------------------------------------------------------



C1 = (SalpTPtable.TP(find(SalpTPtable.Cycle==1)));
C2 = (SalpTPtable.TP(find(SalpTPtable.Cycle==2)));
C4 = (SalpTPtable.TP(find(SalpTPtable.Cycle==4)));
C5 = (SalpTPtable.TP(find(SalpTPtable.Cycle==5)));

C1(find(isnan(C1)))=[];
C2(find(isnan(C2)))=[];
C4(find(isnan(C4)))=[];
C5(find(isnan(C5)))=[];

[h,p,ks2stat] = kstest2(C1,C2)
if h==1
    fn = 'ManuscriptValues.txt'
    fileID = fopen(fn,'a');
    formatSpec = 'Salps: Trophic position of Cycle 1 was statistically different from Cycle 2 according to Kolmogorov-Smirnov: p =  %4.4f  \n';
    fprintf(fileID,formatSpec,[p])
    fclose(fileID)
end
[h,p,ks2stat] = kstest2(C1,C4)
if h==1
    fn = 'ManuscriptValues.txt'
    fileID = fopen(fn,'a');
    formatSpec = 'Salps: Trophic position of Cycle 1 was statistically different from Cycle 4 according to Kolmogorov-Smirnov: p =  %4.4f  \n';
    fprintf(fileID,formatSpec,[p])
    fclose(fileID)
end
[h,p,ks2stat] = kstest2(C2,C4)
if h==1
    fn = 'ManuscriptValues.txt'
    fileID = fopen(fn,'a');
    formatSpec = 'Salps: Trophic position of Cycle 2 was statistically different from Cycle 4 according to Kolmogorov-Smirnov: p =  %4.4f  \n';
    fprintf(fileID,formatSpec,[p])
    fclose(fileID)
end

 fn = 'ManuscriptValues.txt'
 fileID = fopen(fn,'a');
 formatSpec = 'Mean trophic positions were %4.2f +/- %4.2f, %4.2f +/- %4.2f, and %4.2f +/- %4.2f for cycles 1, 2, and 4. \n';
 fprintf(fileID,formatSpec,[nanmean(C1),nanstd(C1),nanmean(C2),nanstd(C2),nanmean(C4),nanstd(C4)])
 fclose(fileID)

ALL = SalpTPtable;
ALL(find(isnan(SalpTPtable.TP)),:)=[];

 fn = 'ManuscriptValues.txt'
 fileID = fopen(fn,'a');
 formatSpec = 'Mean trophic positions for all salps (across cycles) were %4.2f +/- %4.2f for %4.2f to %4.2f mm sized salps. \n';
 fprintf(fileID,formatSpec,[nanmean(SalpTPtable.TP),nanstd(SalpTPtable.TP),nanmin(ALL.Length),nanmax(ALL.Length)])
 fclose(fileID)


C1 = (SizeFracTPtable.TP(find(SizeFracTPtable.Cycle==1)));
C2 = (SizeFracTPtable.TP(find(SizeFracTPtable.Cycle==2)));
C3 = (SizeFracTPtable.TP(find(SizeFracTPtable.Cycle==3)));
C4 = (SizeFracTPtable.TP(find(SizeFracTPtable.Cycle==4)));
C5 = (SizeFracTPtable.TP(find(SizeFracTPtable.Cycle==5)));

C1(find(isnan(C1)))=[];
C2(find(isnan(C2)))=[];
C3(find(isnan(C3)))=[];
C4(find(isnan(C4)))=[];
C5(find(isnan(C5)))=[];
 
[h1,p1,ks2stat] = kstest2(C5,C1)
[h2,p2,ks2stat] = kstest2(C5,C2)
[h3,p3,ks2stat] = kstest2(C5,C3)
[h4,p4,ks2stat] = kstest2(C5,C4)
if h1==1 & h2==1 & h3==1 & h4==1
    fn = 'ManuscriptValues.txt'
    fileID = fopen(fn,'a');
    formatSpec = 'Size Fractionated Zooplankton: Trophic position of Cycle 5 was statistically different from all other cycles according to Kolmogorov-Smirnov: maximum p-value =  %4.4f  \n';
    fprintf(fileID,formatSpec,[max([p1,p2,p3,p4])])
    fclose(fileID)
end

Size1 = (SizeFracTPtable.TP(find(SizeFracTPtable.MedSize>0.2 & SizeFracTPtable.MedSize<0.5)));
Size2 = (SizeFracTPtable.TP(find(SizeFracTPtable.MedSize>0.5 & SizeFracTPtable.MedSize<1)));
Size3 = (SizeFracTPtable.TP(find(SizeFracTPtable.MedSize>1 & SizeFracTPtable.MedSize<2)));
Size4 = (SizeFracTPtable.TP(find(SizeFracTPtable.MedSize>2 & SizeFracTPtable.MedSize<4)));
Size5 = (SizeFracTPtable.TP(find(SizeFracTPtable.MedSize>4)));

Size1(find(isnan(Size1)))=[];
Size2(find(isnan(Size2)))=[];
Size3(find(isnan(Size3)))=[];
Size4(find(isnan(Size4)))=[];
Size5(find(isnan(Size5)))=[];
 
[h1,p1,ks2stat] = kstest2(Size1,Size2)
[h2,p2,ks2stat] = kstest2(Size1,Size3)
[h3,p3,ks2stat] = kstest2(Size1,Size4)
[h4,p4,ks2stat] = kstest2(Size1,Size5)
    fn = 'ManuscriptValues.txt'
    fileID = fopen(fn,'a');
    formatSpec = 'The smallest zooplankton size class had an average trophic position of  %4.4f +/- %4.4f compared to %4.4f +/- %4.4f for all other size classes combined.  \n';
    fprintf(fileID,formatSpec,[mean(Size1),std(Size1),mean([Size2;Size3;Size4;Size5]),std([Size2;Size3;Size4;Size5])])
    fclose(fileID)

    fileID = fopen(fn,'a');
    formatSpec = 'Averaged across all samples, size-fractionated zooplankton had an average trophic position of  %4.4f +/- %4.4f .  \n';
    fprintf(fileID,formatSpec,[mean([Size1;Size2;Size3;Size4;Size5]),std([Size1;Size2;Size3;Size4;Size5])])
    fclose(fileID)
if h1==1
    fn = 'ManuscriptValues.txt'
    fileID = fopen(fn,'a');
    formatSpec = 'Trophic positions were significantly different between size classes 1 and 2 according to Kolmogorov-Smirnov: p =  %4.4f  \n';
    fprintf(fileID,formatSpec,[p1])
    fclose(fileID)
end
if h2==1
    fn = 'ManuscriptValues.txt'
    fileID = fopen(fn,'a');
    formatSpec = 'Trophic positions were significantly different between size classes 1 and 3 according to Kolmogorov-Smirnov: p =  %4.4f  \n';
    fprintf(fileID,formatSpec,[p2])
    fclose(fileID)
end
if h3==1
    fn = 'ManuscriptValues.txt'
    fileID = fopen(fn,'a');
    formatSpec = 'Trophic positions were significantly different between size classes 1 and 4 according to Kolmogorov-Smirnov: p =  %4.4f  \n';
    fprintf(fileID,formatSpec,[p3])
    fclose(fileID)
end
if h4==1
    fn = 'ManuscriptValues.txt'
    fileID = fopen(fn,'a');
    formatSpec = 'Trophic positions were significantly different between size classes 1 and 5 according to Kolmogorov-Smirnov: p =  %4.4f  \n';
    fprintf(fileID,formatSpec,[p4])
    fclose(fileID)
end


for cycle = [1,2,4]
    zoop = (SizeFracTPtable.TP(find(SizeFracTPtable.Cycle==cycle)));
    salp = (SalpTPtable.TP(find(SalpTPtable.Cycle==cycle)));
    salp(find(isnan(salp)))=[];
    zoop(find(isnan(zoop)))=[];
    [h,p,ks2stat] = kstest2(zoop,salp);
    diff = mean(zoop) - mean(salp);
    if h==1
        fileID = fopen(fn,'a');
        formatSpec = 'For cycle %4.1f salp trophic positions were significantly different from size-fractionated trophic positions according to Kolmogorov-Smirnov: p =  %4.4f with an average trophic position difference of %4.4f   \n';
        fprintf(fileID,formatSpec,[cycle,p,diff])
        fclose(fileID)
    else
        fileID = fopen(fn,'a');
        formatSpec = 'For cycle %4.1f salp trophic positions were an average of %4.4f lower than size-fractionated zooplankton, although this difference was not statistically significant according to Kolmogorov-Smirnov: p =  %4.4f  \n';
        fprintf(fileID,formatSpec,[cycle,p,diff])
        fclose(fileID)
    end
end
    

