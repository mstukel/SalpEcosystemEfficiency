function [output] = MakeLegend(Bounds,Shapes,Colors,Labels,FtSize)


ht = Bounds(2,2)-Bounds(2,1);
wid = Bounds(1,2)-Bounds(1,1);
numShapes = height(Labels);
dht = ht/(numShapes*2);
ylocs = [Bounds(2,1)+dht:dht*2:Bounds(2,2)]
xlocs = [Bounds(1,1)+wid*0.02,Bounds(1,1)+wid*0.18,Bounds(1,1)+wid*0.2]

%plot([Bounds(1,1),Bounds(1,2),Bounds(1,2),Bounds(1,1),Bounds(1,1)],[Bounds(2,1),Bounds(2,1),Bounds(2,2),Bounds(2,2),Bounds(2,1)],'-k')  %Outer box
fill([Bounds(1,1),Bounds(1,2),Bounds(1,2),Bounds(1,1),Bounds(1,1)],[Bounds(2,1),Bounds(2,1),Bounds(2,2),Bounds(2,2),Bounds(2,1)],'w')  %Outer box
for i=1:numShapes
    tmp = char(Shapes(i));
    if strcmp(tmp(1),'-')
        plot([xlocs(1),xlocs(2)],[ylocs(i) ylocs(i)],tmp,'Color',Colors(i,:),'LineWidth',2)
    elseif strcmp(tmp(1),'o') | strcmp(tmp(1),'s') | strcmp(tmp(1),'d') | strcmp(tmp(1),'^') 
        plot([xlocs(1)+xlocs(2)]/2,[ylocs(i)],tmp,'Color',Colors(i,:),'MarkerFaceColor',Colors(i,:),'MarkerEdgeColor','k')
    else
        plot([xlocs(1)+xlocs(2)]/2,[ylocs(i)],tmp,'Color',Colors(i,:))
    end
    text(xlocs(3),ylocs(i),Labels(i),'FontSize',FtSize)
end

output = 0;