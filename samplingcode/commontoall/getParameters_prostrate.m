function pars = getParameters_prostrate(planC,caseNum)
% This function prepares all parameters for optimization
% The key of this implementation is for different structures name we assign
% them to different groups and give them different value

initFile = sprintf('init_case_%d.m',caseNum);
fid = fopen(initFile,'w');
stringStepI = 'rectum';
stringStepII = 'rectum';



indexS = planC{end};

structures = planC{indexS.structures};
numStructures = length(structures);

T = [];
RI = [];
RII = [];
RIII = [];
Dpre = [];


structures = changeStructuresName(structures);
k=1;
sampleRates=[];
choppedIM = size(planC{planC{end}.IM},2);
IM = planC{planC{end}.IM}(choppedIM).IMDosimetry;
beamlets = [IM.beams(:).beamlets];

for i = 1:numStructures
    structureName = structures(i).structureName;
    structureName
    if strcmp(lower(structureName),stringStepII)
        fprintf(fid,'%s=%d;',structureName,i);
        fprintf(fid,'\n');
        sampleRates(i) = beamlets(i).sampleRate;
        RII = [RII i];
    end
    
    if strcmp(lower(structureName),stringStepI)
        fprintf(fid,'%s=%d;',structureName,i);
        fprintf(fid,'\n');
        RI = [RI i];
        sampleRates(i) = beamlets(i).sampleRate;
    elseif ~isempty(strfind(structureName,'CTV')) || ...
            (~isempty(strfind(upper(structureName),'PTV')) && (strfind(upper(structureName),'PTV')==1))
        
        fprintf(fid,'PTV%d=%d;',k,i);
        fprintf(fid,'\n');
        
        k=k+1;
        T = [T i];
        dosenum = size(planC{planC{end}.dose},2);
        DEFAULTPRESDOES = Dx(planC, i, dosenum, 80);
        if ~isempty(strfind(structureName,' '))
            structureName(~ismember(structureName, '0':'9')) = '';
            presDose = str2num(structureName)/100;
            Dpre = [Dpre presDose];
        else
            Dpre = [Dpre DEFAULTPRESDOES];
        end
        sampleRates(i) = beamlets(i).sampleRate;
    elseif ~strcmp(structureName,'OUT')
        fprintf(fid,'%s=%d;',structureName,i);
        fprintf(fid,'\n');
        RIII = [RIII i];
        sampleRates(i) = beamlets(i).sampleRate;
    end
end
DMax = max(Dpre)*1.15;
weightsPrescription=ones(length(T));
writeList(fid,'targets',T);
writeList(fid,'step1OARs',RI);
writeList(fid,'step2OARs',RII);
writeList(fid,'step3OARs',RIII);
writeList(fid,'sampleRates',sampleRates);
writeList(fid,'weightsPrescription',weightsPrescription);
writeFloatList(fid,'Dpre',Dpre);
fprintf(fid,'DMax=%f;\n',DMax);
fprintf(fid,'alpha={[0.01]};\n');
fprintf(fid,'step2XWeights={[1]}');
fprintf(fid,'\n');
fprintf(fid,'upperBound=42.7559;');
fprintf(fid,'\n');


% writeTarget(fid,'targets=[',T);
% fprintf(fid,'weightsPrescription=[');
% for i=1:length(T)
%     if(length(T)==1)
%         k=i;
%     elseif(i==length(T))
%         k=1;
%     else
%         k=i+1;
%     end
%     fprintf(fid,'%d',Dpre(i)/Dpre(k));
%     if(i~=length(T))
%         fprintf(fid,',');
%     else
%         fprintf(fid,'];');
%         fprintf(fid,'\n');
%     end
% end
% writePrescribedDose(fid,T,Dpre);
% 
% 
% writeOrgans(fid,'step1OARs=[',RII,structures);
% fprintf(fid,'\n');
% fprintf(fid,'DMax=[%d];',max(Dpre)*1.15);
% fprintf(fid,'\n');
% writeOrgans(fid,'step2OARs=[',RII,structures);
% fprintf(fid,'\n');
% writeOrgans(fid,'step3OARs=[',RIII,structures);
% fprintf(fid,'\n');
% fclose(fid);
% commnd = sprintf('svn add %s',initFile);
% system(commnd);
% commnd = sprintf('svn commit %s -m added',initFile);
% system(commnd);
% 
% 





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
function writeOrgans(fid,str,o,structures)
fprintf(fid,str);
for i=1:length(o)
    structureName = structures(o(i)).structureName;
    str = sprintf('caseParams.%s',structureName);
    fprintf(fid,str);
    if(i~=length(o))
        fprintf(fid,',');
    else
        fprintf(fid,'];');
        fprintf(fid,'\n');
    end
end
end
function writePrescribedDose(fconfig,T,Dpre)
fprintf(fconfig,'dosePrescribedTarget=[');
for i=1:length(T)
    
    fprintf(fconfig,'%d',Dpre(i));
    if(i~=length(T))
        fprintf(fconfig,',');
    else
        fprintf(fconfig,'];');
        fprintf(fconfig,'\n');
    end
end
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
function writeTarget(fid,str,T)
fprintf(fid,str);
for i=1:length(T)
    str = sprintf('caseParams.PTV%d',i);
    fprintf(fid,str);
    if(i~=length(T))
        fprintf(fid,',');
    else
        fprintf(fid,'];');
        fprintf(fid,'\n');
    end
end
end

function structures = changeStructuresName(structures)

for i = 1:length(structures)
    structureName = structures(i).structureName;
    if ~isempty(strfind(structureName,'gtv')) || ...
            ~isempty(strfind(structureName,'GTV')) || ...
            ~isempty(strfind(structureName,'PV')) || ...
            ~isempty(strfind(structureName,'Avoidance')) || ...
            ~isempty(strfind(structureName,'-')) || ...
            ~isempty(strfind(structureName,'_')) || ...
            ~isempty(strfind(structureName,'N3')) || ...
            (i>=2 && strcmp(structureName,structures(i-1).structureName)) || ...
            (~isempty(strfind(upper(structureName),'PTV')) && (strfind(upper(structureName),'PTV')~=1))
        structures(i).structureName = 'OUT';
    else
        structures(i).structureName =strrep(structureName,' ','_');
        if(strcmp(structures(i).structureName,'New_Structure')~=0)
            structures(i).structureName = 'OUT';
        end
    end
    
end
end

