function pos =  getIntensityDiffLeft(x,y,j)
    yindex =  find(y==y(j));
    xatyindex = x(yindex);
    righx = max(xatyindex>x(j));
    for k=1:beamletNum
       if x(k)==righx && y(k)== y(j)
           if pos ~= j
                 pos = k;
           end
       end
    end
end
