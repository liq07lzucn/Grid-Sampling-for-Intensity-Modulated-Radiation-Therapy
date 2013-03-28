function [ bin ] = getpVoxelsGrid( indices,x,y,z,oar,p)
%Get the boundar indicies
%   Detailed explanation goes here

global planC;
bin=[];
index = unique(indices(:,3));
minelement = p;
MAXGRIDSIZE = 100;
MINGRIDSIZE = 1;
remain=0;
for i=1:size(index)
    sliceindex = [];
    %Gives 1 for the slice(index(i) and 0 in other place.
    in = indices(:,3)==index(i);
    %Get indicies corresponding to the slice index(i)
    in = reshape(indices([in in in]), [], 3);
    %Get yindex correspondigetng to slice index(i).
    yval = y(in(:,2));
    xval = x(in(:,1));
    gridSize =round(size(yval,2)/minelement);
    
    if gridSize>MAXGRIDSIZE
        gridSize = MAXGRIDSIZE;
    end
    if gridSize <= MINGRIDSIZE
        row = 1;
        col =1;
        gridSize=1;
    else
        [row,col] = getGridSize(xval,yval,gridSize);
    end
    [xdiff,ydiff] = getXYRange(xval,yval);
    intervalx = xdiff/col;
    intervaly = ydiff/row;
    xmin = min(xval);
    ymin = min(yval);
    xmax=xmin;
    ymax=ymin;
    gridElment = size(yval,2)/gridSize;
    
    t = gridElment * p/100;
    t = round(t);
    %ymxD = max(max(yval),max(xval));
    xmaxD = max(xval);
    ymaxD = max(yval);
    %xmaxD = max(max(yval),max(xval));
   
    for j=1:row
        ymin = ymax;
        ymax = ymin+intervaly;
        xmin = min(xval);
        ymin = floor(ymin);
        xmin = floor(xmin);
        hold on;
        
        for k=1:col
            
           if(j==1)
            if(xmax<xmaxD)   
                line([xmax,xmax],[ymin,ymaxD],'LineWidth',2);
            end
           end
            xmax = xmin+intervalx;
            xmax = ceil(xmax);
            ymax = ceil(ymax);
           
            %plot([xmax,xmax],[ymin,ymxD],'-s');
           
            
            [xpop,ypop] = getPopulation(xval,yval,xmin,ymin,xmax,ymax);
            elementToDraw=t;
            if size(xpop,2) > 0
                if t> size(xpop,2)
                    for q=1:size(xpop,2)
                        [xind,yind] = getIndex(xpop(q),ypop(q),x,y);
                        sliceindex = [sliceindex;[xind yind index(i)]];
                    end
                    remain = remain + (t-size(xpop,2));
                else
                    canadd = (size(xpop,2)-t);
                    if remain>0 && canadd>remain
                        elementToDraw = t + remain;
                        remain = 0;
                    elseif remain>0
                        elementToDraw = t + canadd;
                        remain =  remain-canadd;
                    end
                    
                    randIndex = getUniqueRandValues(size(xpop,2),elementToDraw);
                    for q=1:elementToDraw
                        m = randIndex(q);
                        [xind,yind] = getIndex(xpop(m),ypop(m),x,y);
                        sliceindex = [sliceindex;[xind yind index(i)]];
                    end
                end
                
            end
            %plot(x(sliceindex(:,1)), y(sliceindex(:,2)),'>','MarkerSize',10);
            
            xmin = xmax;
        end
        if(j==1)
            %axis([floor(min(xval)),xmax,min(yval),ymxD+0.5]);
            line([xmaxD,xmaxD],[ymin,ymaxD],'LineWidth',2);
        end
        line([min(xval),xmaxD],[ymin,ymin],'LineWidth',2);
        
        
    end
    line([min(xval),xmaxD],[ymaxD,ymaxD],'LineWidth',2);
    
    
%     totElementToDraw = size(yval,2) * p/100;
%     elementToDraw = round(totElementToDraw- length(sliceindex));
%     if(elementToDraw>0)
%      
%         remain = remain + elementToDraw;
%     end
    
%     if(totElementToDraw>length(sliceindex))
%         nsIndex = getRemainingIndices(yval,y,sliceindex(:,2))
%         elementToDraw = round(totElementToDraw- length(sliceindex));
%         randIndex = getUniqueRandValues(length(nsIndex),elementToDraw);
%         for q=1:elementToDraw
%             m = nsIndex(randIndex(q));
%             [xind,yind] = getIndex(xval(m),yval(m),x,y);
%             sliceindex = [sliceindex;[xind yind index(i)]];
%         end
%     end
    
   bin = [bin;sliceindex];
   if(size(sliceindex,1)>0 && size(sliceindex,2)>0)
        xval = x(in(:,1));
        innerx = x(sliceindex(:,1));
        innery = y(sliceindex(:,2));
        %plot(xval,yval,'b.','MarkerSize',30);
%         hold on;
%         plot(innerx,innery,'r.','MarkerSize',30);
%         hold off;
%         axisHandle = xlabel('X');
%         set(axisHandle,'FontSize',14);
%         axisHandle=ylabel('Y');
%         set(axisHandle,'FontSize',14);
%         tstr = sprintf('Total # of voxels:%d, # of selected voxels:%d\n Sampling Rate:%d\n',size(xval,2),size(sliceindex,1),p);
        %title_handle=title(tstr);
        %set(title_handle,'FontSize',16);
        %set(gca,'FontSize',14);
       
%         dir = '/Users/parastiwari/research/trnk/graphs';
%         if(exist(dir) ~=0)
%             fname = sprintf('%s/fig%d.pdf',dir,p);
%             saveas(gcf, fname, 'jpg');
%         end
        %saveas(gcf, fname, 'png');
   end
end
end
% function [notSelectedIndex]=getRemainingIndices(yval,y,selectedIndex)
%     newy = y(1:size(yval,2));
%     totalIndex = find(ismember(y,yval));
%     notSelectedIndex = setdiff(totalIndex,selectedIndex);
% end
function [val]=getUniqueRandValues(n,m)

r=randperm(n);
val=r(1:m);

end
function [row,col]=getGridSize(x,y,num)

[xdiff,ydiff] = getXYRange(x,y);

if(xdiff>ydiff)
    time = xdiff/ydiff;
    row = round(sqrt(num/time));
    if row==0
        row =1;
    end
    
    col = round(num/row);
    %         if siz>=2
    %             col = floor(siz-row);
    %         else
    %             col = floor(num/row);
    %         end
else
    time = ydiff/xdiff;
    col = round(sqrt(num/time));
    if col ==0
        col =1;
    end
    row = round(num/col);
    %         siz = sqrt(num);
    %         if siz>=2
    %             row = floor(siz-col);
    %         else
    %             row = floor(num/col);
    %         end
end
end
function [xdiff,ydiff]=getXYRange(xval,yval)
xmin = min(xval);
xmax = max(xval);
xdiff = xmax-xmin;
ymin = min(yval);
ymax = max(yval);
ydiff = ymax-ymin;
end
function [xpop,ypop]=getPopulation(xval,yval,xmin,ymin,xmax,ymax)

j =  xval >= xmin & xval<xmax;
k = yval>=ymin & yval <ymax ;
m = j&k;
xpop = xval(m);
ypop = yval(m);
end

function [xind,yind] = getIndex(xval,yval,x,y)
xind = find(x==xval);
yind = find(y==yval);

end
function drawGrid(xval,yval,row,col)
shift =0;
gx = min(xval)-shift;
gy = min(yval)-shift;
[xdiff,ydiff] = getXYRange(xval,yval);
intervalx = xdiff/col;
intervaly = ydiff/row;
for j=1:col+1
    gxv = linspace(gx,gx,100);
    gyv = linspace(gy,max(yval)+shift,100);
    plot(gxv,gyv);
    gx = gx+intervalx;
end


gx = min(xval)-shift;
gy = min(yval)-shift;
for j=1:row+1
    gxv = linspace(gx,max(xval)+shift,100);
    gyv = linspace(gy,gy,100);
    plot(gxv,gyv);
    gy = gy+intervaly;
end
end