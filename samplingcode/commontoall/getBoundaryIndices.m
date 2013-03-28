function [ bin,rows ] = getBoundaryIndices( indices,x,y,z,oar)
%Get the boundar indicies 
%   Detailed explanation goes here

    global planC;
%    dir = 'c:\\IMRT\\figure\\';
    
 %   dir = sprintf('%s%d',dir,oar);
%     stat = rmdir(dir,'s');
%     mkdir(dir);
   
	minslice = min(indices(:,3));
    maxslice = max(indices(:,3));
    index = indices(:,3)==minslice;
    %Get the indicies with the minimum z
    bin = reshape(indices([index index index]), [], 3);
    index = indices(:,3)==maxslice;
    offset = size(bin,1);
    %Get the indicies with maximum z
    bin = [bin;reshape(indices([index index index]), [], 3)];
    index = unique(indices(:,3));
    epsilon = 1;
    oind=[];
   
    for i=2:size(index)-1
	     sliceindex = [];
        %Gives 1 for the slice(index(i) and 0 in other place.
        in = indices(:,3)==index(i); 
        %Get indicies corresponding to the slice index(i)
        in = reshape(indices([in in in]), [], 3);
        %Get yindex corresponding to slice index(i).
        yval = y(in(:,2));
       
         process = zeros(size(y));
        %Iterate for each y
        for j=1:size(yval,2)
            %Flag so that we do not add duplicate element of the same
            %slice. 
           
            yindex = in(j,2);
            
           %Make sure that we haven't already processed this slice.
           if process(yindex) ~= 1 
               process(yindex) = 1;
               xindex = in((in(:,2)==yindex),1);
               %Get all x's at yj.
               xval = x(xindex);
               for k =1:size(xval,2)
                   xind = find(x==xval(k));
                   yinatx = in((in(:,1)== xind),2);
                   yatx = y(yinatx);
                   yvalue = y(yindex);
                   if (lookDown(yvalue-epsilon,yatx(yatx<yvalue)) == 0 ||  lookUp(yvalue+epsilon,yatx(yatx>yvalue))==0|| lookDown(xval(k)-epsilon,xval(xval<xval(k)))==0 || lookUp(xval(k)+epsilon,xval(xval>xval(k)))==0)
                        sliceindex = [sliceindex;[xind yindex index(i)]];
                       %Row number of otherindicies
                        oind =[oind;offset+j+k-1];
                       
                   end
                   %sliceindex = [sliceindex;[xind yindex index(i)]];
%                    if yvalue == max(yatx) || yvalue == min(yatx) || xval(k)==min(xval) || xval(k) == max(xval)
%                         sliceindex = [sliceindex;[xind yindex index(i)]];
%                    end
                   clear xind yinatx yatx yvalue;
               end
               % plot(x(in(:,1)),yval,'d',x(sliceindex(:,1)),y(sliceindex(:,2)),'>'); 
               clear xval  xindex yindex
           end
           
           
        end
        offset = offset+size(yval,2); 
        %if(size(sliceindex,1)>0 && size(sliceindex,2)>0
%         xval = x(in(:,1));
%         bndryx = x(sliceindex(:,1));
%         bndryy = y(sliceindex(:,2));
%         plot(xval,yval,'d',bndryx,bndryy,'>'); 
%         fname = sprintf('%s\\fig%d.png',dir,index(i));
%         saveas(gcf, fname, 'png');
        bin = [bin;sliceindex];
       
        clear xval yval bndryx bndryy sliceindex;
    end
    rows =oind;
end
function indi = getIndicies(val,indicies,in)
      bndryindex = indicies==val;
      indi = reshape(in([bndryindex bndryindex bndryindex]),[],3);
end
function flag =  lookDown(val,arr)
   found = 0;
    for i=1:size(arr,2)
        if arr(i)>=val
            found = 1;
            break;
        end
    end
   flag = found;
    
end
function flag = lookUp(val,arr)
 found = 0;
    for i=1:size(arr,2)
        if arr(i)<=val
            found = 1;
            break;
        end
    end
   flag = found;
end
