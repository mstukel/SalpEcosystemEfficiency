clearvars
close all
warning('off','MATLAB:table:ModifiedAndSavedVarnames')

load('BodyGut.mat')


%---------------------------------------------------------------------------------------------------------------------------------
%---------------Body MINUS Gut Figure------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------


fighandle = figure(1);
fighandle.Units = 'inches';
fighandle.Position = [1 3 3.5 4];


hold on
Labels = {'Bulk','Alanine','Glutamic Acid','Leucine','Aspartic Acid','Isoleucine','Proline','Valine','Glycine','Lysine','Phenylalanine','Serine','Threonine'}
plot([0.5,13.5],[0 0],'-k')
plot([1.5 1.5],[-5 8],'-k','LineWidth',2)
plot([8.5 8.5],[-5 8],'-k','LineWidth',2)
boxplot([BodyGutDiff.x_15NAir___,BodyGutDiff.Ala,BodyGutDiff.Glx,BodyGutDiff.Leu,BodyGutDiff.Asx,BodyGutDiff.Ile,...
       BodyGutDiff.Pro,BodyGutDiff.Val,BodyGutDiff.Gly,BodyGutDiff.Lys,BodyGutDiff.Phe,BodyGutDiff.Ser,BodyGutDiff.Thr],...
       Labels,'PlotStyle','traditional')
ylim([-5 8])
xlim([0.5,13.5])
ylabel('\delta^1^5N Salp body minus gut')
title('         "Trophic" AAs        "Source" AAs')


fn = 'GutPlots.BoxPlot'
exportgraphics(gcf,[fn,'.pdf'],'Resolution',600)
exportgraphics(gcf,[fn,'.png'],'Resolution',600)

TEF_all = [BodyGutDiff.x_15NAir___,BodyGutDiff.Ala,BodyGutDiff.Glx,BodyGutDiff.Leu,BodyGutDiff.Asx,BodyGutDiff.Ile,...
       BodyGutDiff.Pro,BodyGutDiff.Val];
for i=1:width(TEF_all);
    tmp = TEF_all(:,i);
    tmp(find(isnan(tmp)))=[];
    TEF(i,1) = mean(tmp);
    TEF(i,2) = std(tmp);
    TEF(i,3) = std(tmp)./sqrt(length(tmp));
end
TEF2 = array2table(TEF','VariableNames',Labels(1:8),'RowNames',{'Mean','StDev','StErr'});
TEF = array2table(TEF,'VariableNames',{'Mean','StDev','StErr'},'RowNames',Labels(1:8)');
origin = 'GutPlots.m'



%---------------------------------------------------------------------------------------------------------------------------------
%---------------End Body MINUS Gut Figure------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------


%---------------------------------------------------------------------------------------------------------------------------------
%---------------Body v Gut Figure------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------



BodyGut = NaN(290,2);
for i=1:height(Salps)-1
    if Salps.Body1Gut0Hyp2(i)==1
        if Salps.x_15NAir___(i)-Salps.x_15NAir___(i+1)==BodyGutDiff.x_15NAir___(i)
            BodyGut(i,1) = Salps.x_15NAir___(i);
            BodyGut(i,2) = Salps.x_15NAir___(i+1);
        elseif Salps.x_15NAir___(i)-Salps.x_15NAir___(i-1)==BodyGutDiff.x_15NAir___(i)
            BodyGut(i,1) = Salps.x_15NAir___(i);
            BodyGut(i,2) = Salps.x_15NAir___(i-1);
        end
    else
        BodyGut(i,1) = NaN;
        BodyGut(i,2) = NaN;
    end
end
Ala = str2double(Salps.Ala);
Glx = str2double(Salps.Glx);
Leu = str2double(Salps.Leu);
Asx = str2double(Salps.Asx);
Ile = str2double(Salps.Ile);
Pro = str2double(Salps.Pro);
Val = str2double(Salps.Val);
BodyGutAA = NaN(290,14);
for i=2:height(Salps)-1
    if Salps.Body1Gut0Hyp2(i)==1
        if Ala(i)-Ala(i+1)==BodyGutDiff.Ala(i)
            BodyGutAA(i,1) = Ala(i);
            BodyGutAA(i,2) = Ala(i+1);
            BodyGutAA(i,3) = Glx(i);
            BodyGutAA(i,4) = Glx(i+1);
            BodyGutAA(i,5) = Leu(i);
            BodyGutAA(i,6) = Leu(i+1);
            BodyGutAA(i,7) = Asx(i);
            BodyGutAA(i,8) = Asx(i+1);
            BodyGutAA(i,9) = Ile(i);
            BodyGutAA(i,10) = Ile(i+1);
            BodyGutAA(i,11) = Pro(i);
            BodyGutAA(i,12) = Pro(i+1);
            BodyGutAA(i,13) = Val(i);
            BodyGutAA(i,14) = Val(i+1);
        elseif Ala(i)-Ala(i-1)==BodyGutDiff.Ala(i)
            BodyGutAA(i,1) = Ala(i);
            BodyGutAA(i,2) = Ala(i-1);
            BodyGutAA(i,3) = Glx(i);
            BodyGutAA(i,4) = Glx(i-1);
            BodyGutAA(i,5) = Leu(i);
            BodyGutAA(i,6) = Leu(i-1);
            BodyGutAA(i,7) = Asx(i);
            BodyGutAA(i,8) = Asx(i-1);
            BodyGutAA(i,9) = Ile(i);
            BodyGutAA(i,10) = Ile(i-1);
            BodyGutAA(i,11) = Pro(i);
            BodyGutAA(i,12) = Pro(i-1);
            BodyGutAA(i,13) = Val(i);
            BodyGutAA(i,14) = Val(i-1);
        end
    else
        BodyGutAA(i,1) = NaN;
        BodyGutAA(i,2) = NaN;
    end
end
Cycle = Salps.Cycle;

load('CycleColors.mat')

fighandle = figure(2);
fighandle.Units = 'inches';
fighandle.Position = [1 3 7 7];

subplot(3,3,1)
hold on
plot([0 10],[0 10],'-k')
for i=1:5
    inds = find(Cycle==i);
    plot(BodyGut(inds,2),BodyGut(inds,1),'dk','MarkerFaceColor',cols(i,:))
end
temp = [BodyGut(:,2),BodyGut(:,1)];
temp(find(isnan(temp(:,1))),:)=[];
b = regress(temp(:,2),[temp(:,1),ones(size(temp(:,1)))]);
plot( [min(temp(:,1)) max(temp(:,2))] , [min(temp(:,1)) max(temp(:,2))]*b(1)+b(2), '-r','LineWidth',2,'Color','0.7 0 0.7')
xlabel('Gut \delta^1^5N')
ylabel('Body')
set(gca,'box','on')
title('Bulk')
set(gca,'FontSize',9)
text(9,1,'a','FontSize',10)

subplot(3,3,2)
hold on
plot([0 20],[0 20],'-k')
ind = 1;
for i=1:5
    inds = find(Cycle==i);
    plot(BodyGutAA(inds,ind+1),BodyGutAA(inds,ind),'dk','MarkerFaceColor',cols(i,:))
end
temp = [BodyGutAA(:,ind+1),BodyGutAA(:,ind)];
temp(find(isnan(temp(:,1))),:)=[];
b = regress(temp(:,2),[temp(:,1),ones(size(temp(:,1)))]);
plot( [min(temp(:,1)) max(temp(:,2))] , [min(temp(:,1)) max(temp(:,2))]*b(1)+b(2), '-r','LineWidth',2,'Color','0.7 0 0.7')
xlabel('Gut \delta^1^5N')
ylabel('Body')
set(gca,'box','on')
title('Alanine')
set(gca,'FontSize',9)
text(18,2,'b','FontSize',10)

subplot(3,3,3)
hold on
plot([0 20],[0 20],'-k')
ind = 3;
for i=1:5
    inds = find(Cycle==i);
    plot(BodyGutAA(inds,ind+1),BodyGutAA(inds,ind),'dk','MarkerFaceColor',cols(i,:))
end
temp = [BodyGutAA(:,ind+1),BodyGutAA(:,ind)];
temp(find(isnan(temp(:,1))),:)=[];
b = regress(temp(:,2),[temp(:,1),ones(size(temp(:,1)))]);
plot( [min(temp(:,1)) max(temp(:,2))] , [min(temp(:,1)) max(temp(:,2))]*b(1)+b(2), '-r','LineWidth',2,'Color','0.7 0 0.7')
xlabel('Gut \delta^1^5N')
ylabel('Body')
set(gca,'box','on')
title('Glutamic Acid')
set(gca,'FontSize',9)
text(18,2,'c','FontSize',10)

subplot(3,3,4)
hold on
plot([0 20],[0 20],'-k')
ind = 5;
for i=1:5
    inds = find(Cycle==i);
    plot(BodyGutAA(inds,ind+1),BodyGutAA(inds,ind),'dk','MarkerFaceColor',cols(i,:))
end
temp = [BodyGutAA(:,ind+1),BodyGutAA(:,ind)];
temp(find(isnan(temp(:,1))),:)=[];
b = regress(temp(:,2),[temp(:,1),ones(size(temp(:,1)))]);
plot( [min(temp(:,1)) max(temp(:,2))] , [min(temp(:,1)) max(temp(:,2))]*b(1)+b(2), '-r','LineWidth',2,'Color','0.7 0 0.7')
xlabel('Gut \delta^1^5N')
ylabel('Body')
set(gca,'box','on')
title('Leucine')
set(gca,'FontSize',9)
text(18,2,'d','FontSize',10)

subplot(3,3,5)
hold on
plot([0 20],[0 20],'-k')
ind = 7;
for i=1:5
    inds = find(Cycle==i);
    plot(BodyGutAA(inds,ind+1),BodyGutAA(inds,ind),'dk','MarkerFaceColor',cols(i,:))
end
temp = [BodyGutAA(:,ind+1),BodyGutAA(:,ind)];
temp(find(isnan(temp(:,1))),:)=[];
b = regress(temp(:,2),[temp(:,1),ones(size(temp(:,1)))]);
plot( [min(temp(:,1)) max(temp(:,2))] , [min(temp(:,1)) max(temp(:,2))]*b(1)+b(2), '-r','LineWidth',2,'Color','0.7 0 0.7')
xlabel('Gut \delta^1^5N')
ylabel('Body')
set(gca,'box','on')
title('Aspartic Acid')
set(gca,'FontSize',9)
text(18,2,'e','FontSize',10)

subplot(3,3,6)
hold on
plot([0 20],[0 20],'-k')
ind = 9;
for i=1:5
    inds = find(Cycle==i);
    plot(BodyGutAA(inds,ind+1),BodyGutAA(inds,ind),'dk','MarkerFaceColor',cols(i,:))
end
temp = [BodyGutAA(:,ind+1),BodyGutAA(:,ind)];
temp(find(isnan(temp(:,1))),:)=[];
b = regress(temp(:,2),[temp(:,1),ones(size(temp(:,1)))]);
plot( [min(temp(:,1)) max(temp(:,2))] , [min(temp(:,1)) max(temp(:,2))]*b(1)+b(2), '-r','LineWidth',2,'Color','0.7 0 0.7')
xlabel('Gut \delta^1^5N')
ylabel('Body')
set(gca,'box','on')
title('Isoleucine')
set(gca,'FontSize',9)
text(18,2,'f','FontSize',10)

subplot(3,3,7)
hold on
plot([0 20],[0 20],'-k')
ind = 11;
for i=1:5
    inds = find(Cycle==i);
    plot(BodyGutAA(inds,ind+1),BodyGutAA(inds,ind),'dk','MarkerFaceColor',cols(i,:))
end
temp = [BodyGutAA(:,ind+1),BodyGutAA(:,ind)];
temp(find(isnan(temp(:,1))),:)=[];
b = regress(temp(:,2),[temp(:,1),ones(size(temp(:,1)))]);
plot( [min(temp(:,1)) max(temp(:,2))] , [min(temp(:,1)) max(temp(:,2))]*b(1)+b(2), '-r','LineWidth',2,'Color','0.7 0 0.7')
xlabel('Gut \delta^1^5N')
ylabel('Body')
set(gca,'box','on')
title('Proline')
set(gca,'FontSize',9)
text(18,2,'g','FontSize',10)

subplot(3,3,8)
hold on
plot([0 20],[0 20],'-k')
ind = 13;
for i=1:5
    inds = find(Cycle==i);
    plot(BodyGutAA(inds,ind+1),BodyGutAA(inds,ind),'dk','MarkerFaceColor',cols(i,:))
end
temp = [BodyGutAA(:,ind+1),BodyGutAA(:,ind)];
temp(find(isnan(temp(:,1))),:)=[];
b = regress(temp(:,2),[temp(:,1),ones(size(temp(:,1)))]);
plot( [min(temp(:,1)) max(temp(:,2))] , [min(temp(:,1)) max(temp(:,2))]*b(1)+b(2), '-r','LineWidth',2,'Color','0.7 0 0.7')
xlabel('Gut \delta^1^5N')
ylabel('Body')
set(gca,'box','on')
title('Valine')
set(gca,'FontSize',9)
text(18,2,'h','FontSize',10)

subplot(3,3,9)
hold on
plot(0.1, 0.9, 'dk', 'MarkerFaceColor', cols(1,:))
text(0.2, 0.9, 'C1 - SA-Sc', 'FontSize',9)
plot(0.1, 0.8, 'dk', 'MarkerFaceColor', cols(2,:))
text(0.2, 0.8, 'C2 - SA', 'FontSize',9)
plot(0.1, 0.7, 'dk', 'MarkerFaceColor', cols(3,:))
text(0.2, 0.7, 'C3 - ST', 'FontSize',9)
plot(0.1, 0.6, 'dk', 'MarkerFaceColor', cols(4,:))
text(0.2, 0.6, 'C4 - ST', 'FontSize',9)
plot(0.1, 0.5, 'dk', 'MarkerFaceColor', cols(5,:))
text(0.2, 0.5, 'C5 - SA', 'FontSize',9)
plot([0.02 0.17],[0.4 0.4],'-k')
text(0.2, 0.4, '1:1 line', 'FontSize',9)
plot([0.02 0.17],[0.3 0.3], '-r','LineWidth',2,'Color','0.7 0 0.7')
text(0.2, 0.3, 'Regression Line', 'FontSize',9)
ylim([0.2 1])
xlim([-0.05 1])
set(gca,'box','on')
set(gca,'XTick',[])
set(gca,'YTick',[])

fn = 'GutPlots.BodyvGut'
exportgraphics(gcf,[fn,'.pdf'],'Resolution',600)
exportgraphics(gcf,[fn,'.png'],'Resolution',600)

%---------------------------------------------------------------------------------------------------------------------------------
%---------------End Body v Gut Figure------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------



%---------------------------------------------------------------------------------------------------------------------------------
%---------------TDF Calculations------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------



for i=1:length(SrcAA)
    SRC(:,i) = eval([char('BodyGutDiff.'),char(SrcAA(i))]);
    tmp = SRC(:,i);
    tmp(find(isnan(tmp)))=[];
    SRC_se(i) = std(tmp)/sqrt(length(tmp));
end
Src_mean = mean(nanmean(SRC));
Src_se = sqrt(sum(SRC_se.^2))/length(SRC_se);
for i=1:length(TrAA)
    TR(:,i) = eval([char('BodyGutDiff.'),char(TrAA(i))]);
    tmp = TR(:,i);
    tmp(find(isnan(tmp)))=[];
    TR_se(i) = std(tmp)/sqrt(length(tmp));
end
Tr_mean = mean(nanmean(TR));
Tr_se = sqrt(sum(TR_se.^2))/length(TR_se);
TDF_salp = Tr_mean - Src_mean;
TDF_salp_se = sqrt(Src_se.^2 + Tr_se.^2);

['Calculated Trophic Discrimination Factor for salps was ',num2str(TDF_salp,2),' +/- ',num2str(TDF_salp_se,2)]

fn = 'ManuscriptValues.txt'
fileID = fopen(fn,'a');
%fprintf(fileID,['Calculated Trophic Discrimination Factor for salps was ',num2str(TDF_salp,2),' +/- ',num2str(TDF_salp_se,2)])
formatSpec = 'we calculate a TDF for salps (TDFsalp) of  %4.2f +/- %8.3f \n';
fprintf(fileID,formatSpec,[TDF_salp,TDF_salp_se])
fclose(fileID)



    Ala = BodyGutDiff.Ala;
    Glx = BodyGutDiff.Glx;
    Phe = BodyGutDiff.Phe;
    Ala_Phe = Ala - Phe;
    Glx_Phe = Glx - Phe;
    Ala_Phe(find(isnan(Ala_Phe)))=[];
    Glx_Phe(find(isnan(Glx_Phe)))=[];
    Ala_Phe_salp_mean = mean(Ala_Phe);
    Ala_Phe_salp_se = std(Ala_Phe)/sqrt(length(Ala_Phe));
    Glx_Phe_salp_mean = mean(Glx_Phe);
    Glx_Phe_salp_se = std(Glx_Phe)/sqrt(length(Glx_Phe));


