function [ data ] = getData( cases,filePrefix,rate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ydata=zeros(length(cases),length(rate));

for i=1:length(cases)
    
    for j=1:length(rate)
        filename = sprintf('%s%d_%d.txt',filePrefix,cases(i),rate(j));
        fid = fopen(filename);
        while ~feof(fid)

            curr = fscanf(fid,'%f');
            ydata(i,j) = curr/60;
        end
        fclose(fid);
    end
end
data = ydata;


end

