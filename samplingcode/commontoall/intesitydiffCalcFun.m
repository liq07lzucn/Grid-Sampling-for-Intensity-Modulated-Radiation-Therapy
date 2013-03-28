function pos =  getIntensityDiffUp(x,y,j)
     xindex = find(x==x(j));
       yatxindex = y(xindex);
    yabovecurry = min(yatxindex>y(j)); 
    for k=1:beamletNum
           if x(k)==x(xindex)&& y(k)== yabovecurry
                if pos ~= j
                    pos = k;
                end
           end
    end
end




