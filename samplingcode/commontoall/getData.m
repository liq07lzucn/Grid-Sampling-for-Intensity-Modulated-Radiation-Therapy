function [ data ] = getData( cases,filePrefix )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ydata=zeros(length(cases),1);

for i=1:length(cases)
    filename = sprintf('%s%d.txt',filePrefix,cases(i));
    filename
    fid = fopen(filename);
    while ~feof(fid)
        
        curr = fscanf(fid,'%f');
        ydata(i) = curr/60;
    end
end
data = ydata;
fclose(fid);

end

