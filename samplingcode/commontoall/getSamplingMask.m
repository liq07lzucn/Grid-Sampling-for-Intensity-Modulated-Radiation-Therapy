function [ mask ] = getSamplingMask( strucNum,mask3M,rate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sV  = size(mask3M);

if ~exist('planC')
    global planC;
end
indexS = planC{end};
[scanNum, relStructNum] = getStructureAssociatedScan(strucNum, planC);
[x,y,z] = getUniformScanXYZVals(planC{indexS.scan}(scanNum));

linearIndex = find(mask3M ~= 0);
[r c s] =  ind2sub(size(mask3M), linearIndex);
indices = [r c s];
[boundaryIndices,rows] = getBoundaryIndices(indices,x,y,z,strucNum);
remainingIndices = setdiff(indices,boundaryIndices,'rows');
innerIndicies = getpVoxelsGrid(remainingIndices,x,y,z,strucNum,rate);
%len = length(remainingIndices)*rate/100;
%len = round(len);
%innerIndicies = remainingIndices(1:len,:,:);
combinedIndices = [boundaryIndices;innerIndicies];
indexV = sub2ind(size(mask3M),double(combinedIndices(:,1)),double(combinedIndices(:,2)),double(combinedIndices(:,3)));
mask = repmat(logical(0), sV);
mask(indexV)=1;
end

