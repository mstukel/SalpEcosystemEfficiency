
clearvars
close all
warning('off','MATLAB:table:ModifiedAndSavedVarnames')
load('Bulk.mat')

%---------------------------------------------------------------------------------------------------------------------------------
%---------------Salps & Zoops + POM---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------

fighandle = figure(9);
fighandle.Units = 'inches';
fighandle.Position = [3 1 8.3 4];
subplot(1,3,1)
hold on
for i=1:length(cycles)
    inds = find(SizeFrac.Cycle==i);
    h1(i)=plot(SizeFrac.x_13CVPDB___(inds),SizeFrac.x_15NAir___(inds),'ok','MarkerFaceColor',cols(i,:),'MarkerSize',5)
end
for i=1:length(cycles)
    inds = find(POM.Cycle==i);
    tmp = boundary(POM.C13(inds),POM.N15(inds),0);
    h=fill(POM.C13(inds(tmp)),POM.N15(inds(tmp)),'r','LineStyle','none')
    set(h,'FaceAlpha',0.5)
    set(h,'FaceColor',cols(i,:))
end
for i=1:length(cycles)
    inds = find(SizeFrac.Cycle==i & SizeFrac.MedSize==SizeFrac.MedSize(1));
    plot(SizeFrac.x_13CVPDB___(inds),SizeFrac.x_15NAir___(inds),'ok','MarkerFaceColor',cols(i,:),'MarkerSize',8)
    inds = find(SizeFrac.Cycle==i & SizeFrac.MedSize==SizeFrac.MedSize(2));
    plot(SizeFrac.x_13CVPDB___(inds),SizeFrac.x_15NAir___(inds),'dk','MarkerFaceColor',cols(i,:),'MarkerSize',8)
    inds = find(SizeFrac.Cycle==i & SizeFrac.MedSize==SizeFrac.MedSize(3));
    plot(SizeFrac.x_13CVPDB___(inds),SizeFrac.x_15NAir___(inds),'sk','MarkerFaceColor',cols(i,:),'MarkerSize',8)
    inds = find(SizeFrac.Cycle==i & SizeFrac.MedSize==SizeFrac.MedSize(4));
    plot(SizeFrac.x_13CVPDB___(inds),SizeFrac.x_15NAir___(inds),'vk','MarkerFaceColor',cols(i,:),'MarkerSize',8)
    inds = find(SizeFrac.Cycle==i & SizeFrac.MedSize==SizeFrac.MedSize(5));
    plot(SizeFrac.x_13CVPDB___(inds),SizeFrac.x_15NAir___(inds),'pk','MarkerFaceColor',cols(i,:),'MarkerSize',8)
end
xlabel('\delta^1^3C')
ylabel('\delta^1^5N')
title('Size-frac zooplankton')
set(gca,'box','on')
% h=legend([h1(1), h1(2), h1(3), h1(4), h1(5)],{'C1 - SA-Sc','C2 - SA','C3 - ST','C4 - ST','C5 - SA'},'Location','NorthWest')
% set(h,'AutoUpdate','off')
Bounds = [-29.8, -24; 7, 11.8]
Shapes = {'o';'d';'s';'v';'p'}
Colors = [1 1 1;1 1 1; 1 1 1; 1 1 1; 1 1 1]-1;
Labels = {'0.2 - 0.5 mm';'0.5 - 1 mm';'1 - 2 mm';'2 - 4 mm';'>5 mm'}
[output] = MakeLegend(Bounds,Shapes,Colors,Labels,8)
set(gca,'FontSize',9)
text(-29.5,-5.2,'a','FontSize',10)

% fn = 'BulkPlots.Zoop Isotopes - Cycle'
% exportgraphics(gcf,['..\Figures\BulkPlots\',fn,'.pdf'],'Resolution',600)
% exportgraphics(gcf,['..\Figures\BulkPlots\',fn,'.png'],'Resolution',600)

clearvars -except cycles POM PomTable Salps SizeFrac cols


%Salps - Cycle-------------------
figure(9)
subplot(1,3,2)
Salps = Salps(find(Salps.Body1Gut0Hyp2==1),:);

hold on
for i=1:length(cycles)
    inds = find(Salps.Cycle==i);
    plot(Salps.x_13CVPDB___(inds),Salps.x_15NAir___(inds),'ok','MarkerFaceColor',cols(i,:),'MarkerSize',5)
end
for i=1:length(cycles)
    inds = find(POM.Cycle==i);
    tmp = boundary(POM.C13(inds),POM.N15(inds),0);
    h=fill(POM.C13(inds(tmp)),POM.N15(inds(tmp)),'r','LineStyle','none')
    set(h,'FaceAlpha',0.5)
    set(h,'FaceColor',cols(i,:))
end
for i=1:length(cycles)
    inds = find(Salps.Cycle==i);
    plot(Salps.x_13CVPDB___(inds),Salps.x_15NAir___(inds),'ok','MarkerFaceColor',cols(i,:),'MarkerSize',5)
%     inds = find(Salps.Cycle==i & Salps.MedSize==Salps.MedSize(1));
%     plot(Salps.x_13CVPDB___(inds),Salps.x_15NAir___(inds),'ok','MarkerFaceColor',cols(i,:),'MarkerSize',8)
%     inds = find(Salps.Cycle==i & Salps.MedSize==Salps.MedSize(2));
%     plot(Salps.x_13CVPDB___(inds),Salps.x_15NAir___(inds),'dk','MarkerFaceColor',cols(i,:),'MarkerSize',8)
%     inds = find(Salps.Cycle==i & Salps.MedSize==Salps.MedSize(3));
%     plot(Salps.x_13CVPDB___(inds),Salps.x_15NAir___(inds),'sk','MarkerFaceColor',cols(i,:),'MarkerSize',8)
%     inds = find(Salps.Cycle==i & Salps.MedSize==Salps.MedSize(4));
%     plot(Salps.x_13CVPDB___(inds),Salps.x_15NAir___(inds),'vk','MarkerFaceColor',cols(i,:),'MarkerSize',8)
%     inds = find(Salps.Cycle==i & Salps.MedSize==Salps.MedSize(5));
%     plot(Salps.x_13CVPDB___(inds),Salps.x_15NAir___(inds),'pk','MarkerFaceColor',cols(i,:),'MarkerSize',8)
end
xlabel('\delta^1^3C')
ylabel('\delta^1^5N')
title('Salps')
set(gca,'box','on')
ylim([-6 12])
xlim([-30 -18])
% h=legend({'C1 - SA-Sc','C2 - SA','C3 - ST','C4 - ST'},'Location','NorthWest')
% set(h,'AutoUpdate','off')
% set(gca,'FontSize',8)
Bounds = [-29.8, -24; 7, 11.8]
Shapes = {'o';'o';'o';'o';'o'}
Colors = flipud(cols)
Labels = flipud({'C1 - SA-Sc';'C2 - SA';'C3 - ST';'C4 - ST';'C5 - SA'})
[output] = MakeLegend(Bounds,Shapes,Colors,Labels,8)
text(-29.5,-5.2,'b','FontSize',10)



subplot(2,3,3)
hold on
for i=1:length(cycles)
    temp=table2array(POM(find(POM.Cycle==cycles(i)),:));
    depths=unique(temp(:,5));
    clear temp2 temp3
    for j=1:length(depths)
        temp2=temp(find(temp(:,5)==depths(j)),:);
        temp3(j,1)=depths(j);
        temp3(j,2)=mean(temp2(:,10));
        temp3(j,3)=std(temp2(:,10));
        temp3(j,4)=std(temp2(:,10))/sqrt(length(temp2(:,11)));
    end
    temp3(find(temp3(:,3)==0),:)=[];
    for j=1:height(temp3)
        plot([temp3(j,2)-temp3(j,4),temp3(j,2)+temp3(j,4)],[temp3(j,1),temp3(j,1)]+rand*0.5-0.25,'-k','Color',cols(i,:))
    end
        plot(temp3(:,2),temp3(:,1),'-k','Color',cols(i,:))
        plot(temp3(:,2),temp3(:,1),'ok','MarkerFaceColor',cols(i,:),'MarkerEdgeColor','k')
    temp3 = [i*ones(height(temp3),1),temp3];
    if i==1
        POMC13 = temp3;
    else
        POMC13 = [POMC13;temp3];
    end
end
axis ij
ylabel('Depth (m)')
xlabel('\delta^1^3C')
set(gca,'FontSize',9)
set(gca,'box','on')
Bounds = [-28.9, -25.9; 2, 62]
Shapes = {'-r';'-b';'-g';'-m';'-y'}
Colors = cols;
Labels = {'C1 - SA-Sc';'C2 - SA';'C3 - ST';'C4 - ST';'C5 - SA'}
[output] = MakeLegend(Bounds,Shapes,Colors,Labels,8)
ylim([0 126])
xlim([-29 -21.5])
text(-28.8,115,'c','FontSize',10)

subplot(2,3,6)
hold on
for i=1:length(cycles)
    temp=table2array(POM(find(POM.Cycle==cycles(i)),:));
    depths=unique(temp(:,5));
    clear temp2 temp3
    for j=1:length(depths)
        temp2=temp(find(temp(:,5)==depths(j)),:);
        temp3(j,1)=depths(j);
        temp3(j,2)=mean(temp2(:,11));
        temp3(j,3)=std(temp2(:,11));
        temp3(j,4)=std(temp2(:,11))/sqrt(length(temp2(:,11)));
    end
    temp3(find(temp3(:,3)==0),:)=[];
    for j=1:height(temp3)
        plot([temp3(j,2)-temp3(j,4),temp3(j,2)+temp3(j,4)],[temp3(j,1),temp3(j,1)]+rand*0.5-0.25,'-k','Color',cols(i,:))
    end
    plot(temp3(:,2),temp3(:,1),'-k','Color',cols(i,:))
    plot(temp3(:,2),temp3(:,1),'ok','MarkerFaceColor',cols(i,:),'MarkerEdgeColor','k')
    temp3 = [i*ones(height(temp3),1),temp3];
    if i==1
        POMN15 = temp3;
    else
        POMN15 = [POMN15;temp3];
    end
end
axis ij
ylabel('Depth (m)')
xlabel('\delta^1^5N')
set(gca,'FontSize',9)
set(gca,'box','on')
% Bounds = [5, 9.9; 103, 125]
% Shapes = {'-r';'-b';'-g';'-m';'-y'}
% Colors = cols;
% Labels = {'C1 - SA-Sc';'C2 - SA';'C3 - ST';'C4 - ST';'C5 - SA'}
% [output] = MakeLegend(Bounds,Shapes,Colors,Labels,9)
ylim([0 126])
xlim([-3 7.5])
text(-2.6,115,'d','FontSize',10)



fn = 'BulkPlots.Salps & Zoops + POM'
exportgraphics(gcf,[fn,'.pdf'],'Resolution',600)
exportgraphics(gcf,[fn,'.png'],'Resolution',600)


clearvars -except cols POM PomTable Salps SizeFrac

%---------------------------------------------------------------------------------------------------------------------------------
%---------------End Salps & Zoops + POM------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------

origin = 'BulkPlots.m'
%save('POM CycleAve.mat','PomTable','origin')


MLD = [23.4; 48.4; 31.5; 20.7; 21.3];

SalpEnrichment = NaN(height(Salps),1);
SizeFracEnrichment = NaN(height(SizeFrac),1);
for cycle = 1:5
    ind = find(PomTable.Cycle==cycle & PomTable.Depth<MLD(cycle));
    seston = mean(PomTable.N15(ind));
    ind = find(Salps.Cycle==cycle);
    SalpEnrichment(ind,1) = Salps.x_15NAir___(ind) - seston;
    ind = find(SizeFrac.Cycle==cycle);
    SizeFracEnrichment(ind) = SizeFrac.x_15NAir___(ind) - seston;
end
SizeFracEnrichment(isnan(SizeFracEnrichment))=[];
SalpEnrichment(isnan(SalpEnrichment))=[];



fn = 'ManuscriptValues.txt'
fileID = fopen(fn,'a');
formatSpec = 'Size-fractionated zooplankton samples were enriched by %4.1f +/- %4.1f (mean +/- s.d.) relative to mixed-layer seston.\n';
fprintf(fileID,formatSpec,[mean(SizeFracEnrichment),std(SizeFracEnrichment)])
fclose(fileID)


fn = 'ManuscriptValues.txt'
fileID = fopen(fn,'a');
formatSpec = 'Salps were enriched by %4.1f +/- %4.1f (mean +/- s.d.) relative to mixed-layer seston.\n';
fprintf(fileID,formatSpec,[mean(SalpEnrichment),std(SalpEnrichment)])
fclose(fileID)
