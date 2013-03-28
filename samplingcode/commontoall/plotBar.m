function plotBar(X,Y,fname,lgnd,xlbl,ylbl,title,ytick,s)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
   
      
    bar(Y,s);
    set(gca,'XTickLabel',X);
    set(gcf,'name',title);
    mn = ceil(min(min(Y)));
    mx = floor(max(max(Y)));
   
    %set(gca,'YTickLabel',ytick);
    if(ytick)
        set(gca,'YTick',[mn:2:mx]);
    end
    if(size(lgnd,1) ~= 0)
        h_legend=legend('Location', 'NorthEastOutside', lgnd);  
        set(h_legend,'FontSize',14);

    end
    axisHandle=xlabel(xlbl);
    set(axisHandle,'FontSize',18);
    axisHandle=ylabel(ylbl);
    set(axisHandle,'FontSize',18);
    set(gca,'FontSize',14);
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20]);
    
    saveas(gcf, fname, 'jpg');
    
end

