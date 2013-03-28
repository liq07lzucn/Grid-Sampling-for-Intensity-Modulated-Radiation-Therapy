function [ bin ] = getpVoxel( indices,x,y,z,oar,p)
%Get the boundar indicies 
%   Detailed explanation goes here

    global planC;
    directory = 'c:\bIMRT\inner\';
    
    directory = sprintf('%s%d\\%d',directory,oar,p);
    stat = rmdir(directory,'s');
    mkdir(directory);
   
    bin=[];
    index = unique(indices(:,3));
    
    for i=2:size(index)-1
	     sliceindex = [];
        %Gives 1 for the slice(index(i) and 0 in other place.
        in = indices(:,3)==index(i); 
        %Get indicies corresponding to the slice index(i)
        in = reshape(indices([in in in]), [], 3);
        %Get yindex correspondigetng to slice index(i).
        yval = y(in(:,2));
        xval = x(in(:,1));
        num = size(yval,2)*p/100;
%         xind=0;
%         yind =0;
%         [xind,yind] = getIndex(xval(1),yval(1),x,y);
%         if num ~=0
%             sliceindex = [xind,yind,index(i)];
%         end
%         for m=2:num
%             max=0;
%            
%             for j=1:size(yval,2)
%                 min = 99999;
%                 
%                 for k=1:size(sliceindex,1)
%                    
%                     dist = distance(xval(j),yval(j),x(sliceindex(k,1)),y(sliceindex(k,2)));
%                     if dist<min
%                         min = dist;
%                     end
%                 end
%                 if max<min
%                     max = min;
%                    
%                     selectedIndex = j;
%                    
%                 end
%             end
%              
%              [xind,yind] = getIndex(xval(selectedIndex),yval(selectedIndex),x,y);
%              sliceindex = [sliceindex;[xind yind index(i)]];
%              yval(selectedIndex) = [];
%              xval(selectedIndex) = [];
%              max=0;
%              min = 99999;
%         end
        if num ~=0
            k = size(yval,2)/num;
        else
            k=size(yval,2)+1;
        end
        
            
       j=k;
       while j<=(size(yval,2))
           m = j+k.*rand();
           m=floor(m);
           if m< size(yval,2)            
               xind = find(x==xval(m));
               yind = find(y==yval(m));
               sliceindex = [sliceindex;[xind yind index(i)]];
           end
           j= j+k;
        end
        bin = [bin;sliceindex];
        if size(sliceindex,2)~= 0
            innerx = x(sliceindex(:,1));
            innery = y(sliceindex(:,2));
            yval = y(in(:,2));
            xval = x(in(:,1));
            plot(xval,yval,'d',innerx,innery,'>','MarkerSize',10); 
            fname = sprintf('%s\\fig%d.png',directory,index(i));
            saveas(gcf, fname, 'png');
        end
    end
     
    
end
function [xind,yind] = getIndex(xval,yval,x,y)
    xind = find(x==xval);
    yind = find(y==yval);

end
function dist =  distance(x1,y1,x2,y2)
    val = sqrt((x2-x1)^2+(y2-y1)^2);
    dist = val;
end
