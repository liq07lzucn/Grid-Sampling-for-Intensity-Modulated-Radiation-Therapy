%Call prioritized prescription optimization (priPresOpt) for one case
%  Can use similarly to a script by modifying params below, or use with all
%  default parameters by not entering any arguments, or may pass in inputs as
%  function arguments.
%Created by Vanessa Clark, 10/26/07
 
% copyright (c) 2006-2007, Washington University in St. Louis.
% 
% Permission is granted to use or modify only for non-commercial, 
% non-treatment-decision applications, and further only if this header is 
% not removed from any file. No warranty is expressed or implied for any 
% use whatever: use at your own risk.  Users can request use of CERR for 
% institutional review board-approved protocols.  Commercial users can 
% request a license.  Contact Joe Deasy for more information 
% (radonc.wustl.edu@jdeasy, reversed).
 
 
%% default inputs

function priPresOpt_single(caseNum,samplingRate,loadDirectory,saveDirectory,caseType)
if ~exist('quadslip'),quadslip = 1; end
if ~exist('minslip'), minslip = 0.015; end
if ~exist('partNumStr'),partNumStr = ''; end
 
if ~exist('quadslipStep2'), quadslipStep2 = quadslip; end %DEFAULT
if ~exist('quadslipStep3'), quadslipStep3 = 2*quadslip+quadslip^2; end %DEFAULT
if ~exist('quadslipStep4'), quadslipStep4 = quadslip^3+3*quadslip^2+3*quadslip; end %DEFAULT
if ~exist('quadslipStep5'), quadslipStep5 = -quadslip^4 + 4*quadslip^3 - 6*quadslip^2 + 4*quadslip; end %Changed
 
%if ~exist('saveDirectory'), saveDirectory = [cd '/']; end %DEFAULT: working directory. NOTE: may need to change slash to '/'

if ~exist('saveTag'), saveTag = ''; end
if ~exist('savePlanC'), savePlanC = 1; end %whether or not to save planC (prioptStats will still be saved)
patID = caseNum; 
if ~exist('patID'), patID = caseNum; end
if ~exist('initTag'), initTag = ''; end %This is the tag that goes after the case number in the init_pat file name.  Blank is fine.  %'_small_priopt17b';
if ~exist('extension'), extension = ''; end %This goes after the case number in the plan data file name. % '_wIM_7b_wMoat_wOuterHotspot_CERR3';
if ~exist('optFormulation'), optFormulation = 'orig'; end % 'orig' 'N' 'T' 'W'
 
%% initial setup
 
clear global planC
clear global stateS
global planC % These must be global, in order to be used with certain CERR functions.
global stateS
 
%distinguish this plan in file
%loadDirectory = '/research-projects/tantra/tiwarip/imrt/';
%saveDirectory = '/research-projects/tantra/tiwarip/imrt/';
% loadDirectory ='~/imrt/';
% saveDirectory='~/imrt/';
fid = fopen([saveDirectory 'PriOpt_case' num2str(caseNum) initTag '.txt'],'at+');
fprintf(fid,['\nQuadSlip = ' num2str(quadslip) ', MinSlip = ' num2str(minslip) ':\n']);
fclose(fid);

disp('Loading data...');
eval(['load ' loadDirectory 'case',int2str(patID),extension, '.mat']); %_9b_s1_th01
stateS.CERRFile = [loadDirectory 'case',int2str(patID),extension, '.mat'];

stateS.MLVersion = 7.0; 
 
%% do optimization
disp('Optimizing...');

prioptOptions = createPrioptOptions(initTag, quadslip, minslip, quadslipStep2, quadslipStep3, quadslipStep4, optFormulation, extension, saveTag, quadslipStep5);
prioptOptions.samplingRate = samplingRate;
prioptOptions.caseType=caseType;
[prioptStats] = priPresOpt(prioptOptions, saveDirectory);

 

clear lastDose doseNum V40Rectum MOH84Rectum D90PTV1 D95PTV1 
clear numStructs rectumInd PTV1Ind stepNum numSteps  
 

%% clear data (except prioptStats, prioptOptions, caseParams)
clear global stateS
clear quadslipStep2 quadslipStep3 quadslipStep4
clear saveTag initTag quadslip minslip fid patID partNumStr extension
end 
