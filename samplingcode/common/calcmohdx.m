fName = sprintf(['Result.txt']);
f = fopen(fName,'at+');
fprintf(f,'Resulf of Case %d with boundary and 100% of voxel\n\n',patID);
%fprintf(f,'Resulf of Case %d with boundary and 10% of voxel\n',patID);
fprintf(f,'Total time needed for optimization (sec): %g \n\n', tTotal );
lastDose = size(planC{planC{end}.dose},2);
doseNum = lastDose
%k=1;
%len = 4*length(targets)+ 2*length(step2OARs)+3*length(step3OARs);
%dose = zeros(len,1);
fprintf(f,'Beamlets weight\n');
for j=1:length(wV)
    fprintf(f,'%f  ',wV(j));
end
fprintf(f,'\n');
dose=0;
n=0;

for tar=1:length(targets)
    
    d = influenceM(allVoxelC{targets(tar)},:) * wV / doseScale;
    dose = Dx(planC, targets(tar), doseNum, 95);
    fprintf(f,'D95_%s:%6.2f\n',planC{indexS.structures}(targets(tar)).structureName,dose);
    fprintf('D95_%s:%6.2f\n',planC{indexS.structures}(targets(tar)).structureName,dose);
    dose = max(d);
    fprintf(f,'Max_%s:%6.2f\n',planC{indexS.structures}(targets(tar)).structureName,dose);
    fprintf('Max_%s:%6.2f\n',planC{indexS.structures}(targets(tar)).structureName,dose);
    dose = mean(d);  
    fprintf(f,'mean %s:%6.2f\n',planC{indexS.structures}(targets(tar)).structureName,dose);
    fprintf('mean %s:%6.2f\n',planC{indexS.structures}(targets(tar)).structureName,dose);
    dose = min(d);
    fprintf(f,'Min %s:%6.2f\n',planC{indexS.structures}(targets(tar)).structureName,dose);
    fprintf('Min %s:%6.2f\n',planC{indexS.structures}(targets(tar)).structureName,dose);
    
    
end
for j=1:length(step2OARs)
    d = influenceM(allVoxelC{step2OARs(j)},:) * wV / doseScale;
    dose = MOHx(planC, step2OARs(j), doseNum, 5);   
    fprintf(f,'MOH5 %s:%6.2f\n',planC{indexS.structures}(step2OARs(j)).structureName,dose);
    fprintf('MOH5 %s:%6.2f\n',planC{indexS.structures}(step2OARs(j)).structureName,dose);
    dose= max(d);
    fprintf(f,'Max %s:%6.2f\n',planC{indexS.structures}(step2OARs(j)).structureName,dose);
    fprintf('Max %s:%6.2f\n',planC{indexS.structures}(step2OARs(j)).structureName,dose);
    dose= mean(d);  
    fprintf(f,'Mean %s:%6.2f\n',planC{indexS.structures}(step2OARs(j)).structureName,dose);
    fprintf('Mean %s:%6.2f\n',planC{indexS.structures}(step2OARs(j)).structureName,dose);
    dose= min(d);
    fprintf(f,'Min %s:%6.2f\n',planC{indexS.structures}(step2OARs(j)).structureName,dose);
    fprintf('Min %s:%6.2f\n',planC{indexS.structures}(step2OARs(j)).structureName,dose);    
end

maxdose=[];
for j=1:length(step3OARs)
    d = influenceM(allVoxelC{step3OARs(j)},:) * wV / doseScale;
    dose = max(d);
    fprintf(f,'Max %s:%6.2f\n',planC{indexS.structures}(step3OARs(j)).structureName,dose);
    fprintf('Max %s:%6.2f\n',planC{indexS.structures}(step3OARs(j)).structureName,dose);
    dose= mean(d);
    fprintf(f,'Mean %s:%6.2f\n',planC{indexS.structures}(step3OARs(j)).structureName,dose);
    fprintf('Mean %s:%6.2f\n',planC{indexS.structures}(step3OARs(j)).structureName,dose);
    dose = min(d);
    fprintf(f,'Min %s:%6.2f\n',planC{indexS.structures}(step3OARs(j)).structureName,dose);
    fprintf('Min %s:%6.2f\n',planC{indexS.structures}(step3OARs(j)).structureName,dose);  
end
fprintf(f,'\n');
fprintf(f,'\n\n\n');
fprintf('\n');
fclose(f);



