
% This function prepares all parameters for optimization
% The key of this implementation is for different structures name we assign
% them to different groups and give them different value
function getParametersHN(planC,caseNum)
initFile = sprintf('init_case_%d.m',caseNum);
fid = fopen(initFile,'w');
stringStepI = 'Brainstem';
stringStepII = 'Brainstem';
choppedIM = size(planC{planC{end}.IM},2);
indexS = planC{end};

structures = planC{indexS.structures};
numStructures = length(structures);

targets = [];
step1OARs = [];
step2OARs = [];
step3OARs = [];
Dpre = [];
DEFAULTPRESDOES = 70;
sampleRates=[];
structures = changeStructuresName(structures);
IM = planC{planC{end}.IM}(choppedIM).IMDosimetry;
beamlets = [IM.beams(:).beamlets];

for i = 1:numStructures
    structureName = structures(i).structureName;
    if strcmp(structureName,stringStepII) || strcmp(structureName,'Brianstem')
        step2OARs = [step2OARs i];
        writeStructure(fid,structureName,i);
        sampleRates(i) = beamlets(i).sampleRate;
    end
    if strcmp(structureName,stringStepI) || strcmp(structureName,'Brianstem')
        step1OARs = [step1OARs i];
    elseif ~isempty(strfind(structureName,'CTV'))
        targets = [targets i];
        writeStructure(fid,structureName,i);
        if ~isempty(strfind(structureName,'_'))
            structureName(~ismember(structureName, '0':'9')) = '';
            presDose = str2num(structureName)/100;
            Dpre = [Dpre presDose];
        else
            Dpre = [Dpre DEFAULTPRESDOES];
        end
        sampleRates(i) = beamlets(i).sampleRate;
    elseif ~strcmp(structureName,'OUT')
        step3OARs = [step3OARs i];
        sampleRates(i) = beamlets(i).sampleRate;
        writeStructure(fid,structureName,i);
    end
end
DMax = max(Dpre)*1.15;
weightsPrescription=ones(length(targets));
writeList(fid,'targets',targets);
writeList(fid,'step1OARs',step1OARs);
writeList(fid,'step2OARs',step2OARs);
writeList(fid,'step3OARs',step3OARs);
writeList(fid,'sampleRates',sampleRates);
writeList(fid,'weightsPrescription',weightsPrescription);
writeFloatList(fid,'Dpre',Dpre);
fprintf(fid,'DMax=%f;\n',DMax);
fprintf(fid,'alpha={[0.01]};\n');

%alpha   = {[0.01]}; % 0.01 is the best Brainstem V62 best correlated MOHx
fprintf(fid,'step2XWeights={[1]}');
fprintf(fid,'\n');
fprintf(fid,'upperBound=42.7559;');
fprintf(fid,'\n');
fclose(fid);
pause(5);
end




function writeStructure(fid,name,id)
fprintf(fid,name);
fprintf(fid,'=%d;',id);
fprintf(fid,'\n');
end

function writeFloatList(fid,str,lst)
fprintf(fid,str);
fprintf(fid,'=[');
for i=1:length(lst)
    fprintf(fid,'%f',lst(i));
    if(i<length(lst))
        fprintf(fid,',');
    else
        fprintf(fid,'];');
    end    
end
fprintf(fid,'\n');
end

function writeList(fid,str,lst)
fprintf(fid,str);
fprintf(fid,'=[');
for i=1:length(lst)
    fprintf(fid,'%d',lst(i));
    if(i<length(lst))
        fprintf(fid,',');
    else
        fprintf(fid,'];');
    end
end
fprintf(fid,'\n');
end