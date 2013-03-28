
cases = [9];
loadDirectory = '~/imrt/';
extension='';
global planC;
for i=1:size(cases,2)
    caseNum=cases(i);
    disp('Loading data...');
    eval(['load ' loadDirectory 'case',int2str(caseNum),extension, '.mat']); %_9b_s1_th01
    disp('done loading...');
    [caseParams oMap] = getCaseParams(caseNum);
    eval(['init_pat_xx_', int2str(caseNum)]);
    calcDoseNS
end



