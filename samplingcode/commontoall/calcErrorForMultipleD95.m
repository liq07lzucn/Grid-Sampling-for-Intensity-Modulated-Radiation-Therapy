function calcErrorForMultipleD95()
%cases = [4, 5, 9, 10, 11, 12, 13, 14, 15, 16];
%cases = [5, 25, 104, 108, 113, 123, 208, 30, 79];
cases = [5, 25];

rate=[30];
loaddirectory = '~/imrt/final/';
savedirectory = '~/imrt/final/error/';
mkdir(savedirectory);
fid = fopen('result.txt','w');
dataforplot = zeros(4,length(cases));
lgnd={};
for i=1:length(cases)
    for j=1:length(rate)
        pfname = sprintf('%scase%d_%d_res.mat',loaddirectory,cases(i),rate(j));
        [d95 moh meanNormal] = getMetrics(cases(i),pfname);
        dataforplot(1:4,j) = [d95 moh meanNormal];
        
    end
    pfname = sprintf('%scase%d_res_full.mat',loaddirectory,cases(i));
    [d95 moh meanNormal] = getMetrics(cases(i),pfname);
    dataforplot(1:4,length(rate)+1) = [d95 moh meanNormal];
   
end
for i=1:length(rate)
    lgnd{end+1} = num2str(rate);
end
lgnd{end+1} = 'No Sampling';
fname = sprintf('%serrormultiple.png',savedirectory);
plotBar({'PTV 1(D95)' 'PTV2(D95)' 'BrainStem(Moh 5)' 'Normal Organs(mean Dose)'},dataforplot,fname,lgnd,'structure','Dose(Gy)','Dose Distribution.',false);
end