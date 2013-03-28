function prioptOptions = createPrioptOptions(initTag, slipInQuadTerm, slipMinDoseFactor, slipInQuadTermStep2, slipInQuadTermStep3, slipInQuadTermStep4, optFormulation, extension, saveTag, slipInQuadTermStep5)
%This function returns a default prioptOptions object for use as input to
%priopt.
%
%Created by Vanessa H Clark, 3/13/07, last modified by Vanessa Clark 3/08
%
%Slip values are at the default value used in the LAA paper.
%Note that slipInQuadTerm defines all quad slip terms, while specifying
%later slip term values as parameters overrides their values.
 
% copyright (c) 2006-2007, Washington University in St. Louis.
% 
% Permission is granted to use or modify only for non-commercial, 
% non-treatment-decision applications, and further only if this header is 
% not removed from any file. No warranty is expressed or implied for any 
% use whatever: use at your own risk.  Users can request use of CERR for 
% institutional review board-approved protocols.  Commercial users can 
% request a license.  Contact Joe Deasy for more information 
% (radonc.wustl.edu@jdeasy, reversed).
 
 
prioptOptions = struct('prioptCode', [], 'caseInfo', []);
 
%% Priopt code variables
prioptCode = struct();
 
% prioptCode.step4.minmaxdose = 0;
% prioptCode.step4.minsquareddose = 1;
s=1;
prioptCode.slip.quadSlipStep2 = s;
prioptCode.slip.quadSlipStep3 = 2*s + s^2;
prioptCode.slip.quadSlipStep4 = s^3 + 3*s^2 + 3*s;
prioptCode.slip.minSlipStep4 = 0.015;

prioptCode.maxstep =4;
prioptCode.useMaximizationOfMinDose = 1;
prioptCode.lowerLimit = 0.0001;
prioptCode.csObjType = 'mDx';
prioptOptions.prioptCode = prioptCode;
% %also allow target min dose to slip in step 4:
% if ~exist('slipMinDoseFactor')
%     prioptCode.slip.minSlipStep4 = 0.015;
%     warning('Using default value for slipMinDoseFactor of 0.015');
% else
%     prioptCode.slip.minSlipStep4 = 0;%slipMinDoseFactor;
% end
%  
% if ~exist('optFormulation')
%     prioptCode.optFormulation = 'ORIG';%'N';% %optimization formulation type
%     warning('Using default value for optFormulation of ''orig''');
% else
%     prioptCode.optFormulation = optFormulation;
% end
%  
% switch optFormulation
%     case 'orig'
%         prioptCode.maxstep = 4;  % maximum number of optimization steps
%     case {'N','W','T','P'}
%         prioptCode.maxstep = 5;  % maximum number of optimization steps
%         
%         if ~exist('slipInQuadTermStep5')
%             prioptCode.slip.quadSlipStep5 = -slipInQuadTerm^4 + 4*slipInQuadTerm^3 - 6*slipInQuadTerm^2 + 4*slipInQuadTerm;
%             warning('Using default value for slipInQuadTermStep5, same as value for slipInQuadTerm');
%         else
%             prioptCode.slip.quadSlipStep5 = slipInQuadTermStep5;
%         end
%         if ~exist('slipInQuadTermStep6')
%             prioptCode.slip.quadSlipStep6 = slipInQuadTerm^5 - 5*slipInQuadTerm^4 + 10*slipInQuadTerm^3 - 10*slipInQuadTerm^2 + 5*slipInQuadTerm;
%             warning('Using default value for slipInQuadTermStep6, same as value for slipInQuadTerm');
%         else
%             prioptCode.slip.quadSlipStep6 = slipInQuadTermStep6;
%         end
%     otherwise
%         warning('wrong optFormulation code')
% end
%  
% prioptCode.useMultiplePlots = 1; % for plotDVHs.m: 0 means one DVH with all structures, 1 means separate DVH's for all structures
%  
% prioptCode.lowerLimit = 0.0001; % beamlet with a relative weight below this number will be set to 0 after each step
%  
% prioptCode.csObjType = 'mDx'; %could also be 'meandose'
%  
% prioptCode.useMaximizationOfMinDose = 1; % enable maximization of min dose in PTV's
%  
% % slip values per step: 0 means no slip, 0.2 means 20% in variance, i.e.
% % 10% in standard deviation 
%  
% %Note that setting the value for slipInQuadTerm automatically defines the
% %values for quadslip at each step, while setting the values for any other
% %quadslip term overrides the more generic slipInQuadTerm.  The default
% %quadslip values are based off of the original version of slip, which used
% %powers of 2 and 3 in steps 3 and 4, rather than the current power of 1 in
% %all steps.
% if ~exist('slipInQuadTerm')
%     slipInQuadTerm = 1.0;
%     prioptCode.slip.quadSlipStep2 = 0;%1.0;
%     prioptCode.slip.quadSlipStep3 = 0;%2*1.0 + 1.0^2;
%     prioptCode.slip.quadSlipStep4 = 0;%1.0^3 + 3*1.0^2 + 3*1.0;    
%     warning('Using default value for slipInQuadTerm of 1.0');
% else
%     prioptCode.slip.quadSlipStep2 = 0;%slipInQuadTerm;
%     prioptCode.slip.quadSlipStep3 = 0;%2*slipInQuadTerm + slipInQuadTerm^2;
%     prioptCode.slip.quadSlipStep4 = 0;%slipInQuadTerm^3 + 3*slipInQuadTerm^2 + 3*slipInQuadTerm;
% end
% if ~exist('slipInQuadTermStep2')
%     prioptCode.slip.quadSlipStep2 = slipInQuadTerm;
%     warning('Using default value for slipInQuadTermStep2, same as value for slipInQuadTerm');
% else
%     prioptCode.slip.quadSlipStep2 = slipInQuadTermStep2;
% end
% if ~exist('slipInQuadTermStep3')
%     prioptCode.slip.quadSlipStep3 = 2*slipInQuadTerm + slipInQuadTerm^2;
%     warning('Using default value for slipInQuadTermStep3, same as value for slipInQuadTerm');
% else
%     prioptCode.slip.quadSlipStep3 = slipInQuadTermStep3;
% end
% if ~exist('slipInQuadTermStep4')
%     prioptCode.slip.quadSlipStep4 = slipInQuadTerm^3 + 3*slipInQuadTerm^2 + 3*slipInQuadTerm;
%     warning('Using default value for slipInQuadTermStep4, same as value for slipInQuadTerm');
% else
%     prioptCode.slip.quadSlipStep4 = slipInQuadTermStep4;
% end
%  
% prioptOptions.prioptCode = prioptCode;
%  
% %% Case-specific variables
% caseInfo = struct();
%  
% %caseInfo.patID = 0; %automatically detected
% if ~exist('initTag')
%     warning('Using empty initTag');
%     initTag = '';
% end    
% caseInfo.initTag = initTag;
%  
% if ~exist('saveTag')
%     warning('Using empty saveTag');
%     saveTag = '';
% end    
% caseInfo.saveTag = saveTag;
%  
% if ~exist('extension')
%     warning('Using empty extension');
%     extension = '';
% end    
% caseInfo.extension = extension;
%  
% prioptOptions.caseInfo = caseInfo;
% prioptOptions.caseType='HN';
% 
