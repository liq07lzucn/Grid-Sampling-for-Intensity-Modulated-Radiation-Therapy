function [prioptStats,wV, timesV] = priPresOpt(prioptOptions, saveDirectory)
% function to run prioritized prescription optimization for IMRT
% By Vanessa Clark, 10/25/2007, based on code by Jan Wilkens, 01/30/2006
%
% --For info on how to run this:  see attached documentation,
% "How to run prioritized optimization on a plan.doc"
% --For info on the general procedure:  see attached documentation,
% "PrioritizedPrescriptionOptimizationOverview_Prostate.ppt"
%
% Inputs (that could vary regularly):
%init_pat script file (e.g. init_pat_1) with patient-specific options
%planC, stateS (global variables, not passed as input parameters)
%   Note: it is assumed that the input planC file name (in stateS.CERRFile) contains "caseXX", where XX is the patID
%
% Inputs (that should be static in the final version but for testing may not be):
%prioptOptions (optimization-specific options object created by calling a
%function createPrioptOptions)
%init_pat_xx (script file with prescription-specific options)
%
%
% Outputs (explicit):
%prioptStats (inc. resC), with information about dose at each step of optzn
%caseParams (structure numbers inputted in init_pat_xx)
%wV (beamlet intensity values)
%timesV (times, in seconds, for each step of the optimization)
%
% Outputs (implicit / stored in files):
%planC (global variable, with plan information)
%output to Matlab display
%[nameRun extension initTag '_mosek.txt']
%[nameRun extension initTag '.txt']
%

% copyright (c) 2006-2007, Washington University in St. Louis.
%
% Permission is granted to use or modify only for non-commercial,
% non-treatment-decision applications, and further only if this header is
% not removed from any file. No warranty is expressed or implied for any
% use whatever: use at your own risk.  Users can request use of CERR for
% institutional review board-approved protocols.  Commercial users can
% request a license.  Contact Joe Deasy for more information
% (radonc.wustl.edu@jdeasy, reversed).


% NOTE: DoseScale functionality may be broken, by VHC.

% Modified 2/7/08, VHC, changed 'solution' to 'solutions'
% Modified 2/26/08, VHC, saved in 'solution' and 'solutions' for now
tic
slip = prioptOptions.prioptCode.slip;
if ~exist('saveDirectory')
    saveDirectory = '';
end

%--------------------------------------------------------------------------
% patient specific input:

global planC; % (input parameter now)
indexS = planC{end};
global stateS; % (input parameter now)

rate = prioptOptions.samplingRate;
% auto detect patID, i.e. case number
% it is assumed that the file name contains "caseXX", where XX is the patID
patIDpos = max(strfind(stateS.CERRFile,'case'))+4;
patIDendPos = strfind(stateS.CERRFile,'.')-1;
patID = str2double(stateS.CERRFile(patIDpos:patIDendPos));
if isempty(patID)
    patID = 0;
end
clear patIDpos
if (patID <= 0) %(patID < 50) || (patID > 70)
    disp('patID not recognized!');
    return
end
initFileName = sprintf('init_case_%d',patID);
initwdExt = strcat(initFileName,'.m');
if(strcmp(prioptOptions.caseType,'HN')~=0)
    if(exist(initwdExt,'file')==0)
        getParametersHN(planC,patID);
    end
   
elseif(strcmp(prioptOptions.caseType,'LN')~=0)
    if(exist(initwdExt,'file')==0)
        getParameters_lung(planC,patID);
    end
else
     if(exist(initwdExt,'file')==0)
        getParameters_prostrate(planC,patID);
    end
end
eval(initFileName) 
%--------------------------------------------------------------------------
% setup data:

startstep = 1;

clear prioptStats;

batchMode = 0; % 0 means that DVH's and other figures are shown
try
    sliceCallBack('refresh');
    if ~exist('fluenceMapFigure','var')
        fluenceMapFigure = figure;
    end
    if ~exist('beamWeightFigure','var')
        beamWeightFigure = figure;
    end
    if ~exist('dvhFigure','var')
        dvhFigure = figure;
    end
catch
    batchMode = 1; % no figures in batch mode
end

%--------------------------------------------------------------------------
disp('Preparing data...')

%get index vectors for all relevant structures
voxelC = cell(max([targets step2OARs step3OARs]),1);
numVoxels = zeros(1,max([targets step2OARs step3OARs]));
getBoundaryVoxel
voxelC
getAllVoxel
allVoxelC
for i=1:length(targets)
    %voxelC{targets(i)} = allVoxelC{targets(i)};
    %numVoxels(targets(i)) = numAllVoxels(targets(i));
    if length(voxelC{targets(i)})==0
        fprintf('Target %d has zero voxel',targets(i));
        return;
        
    end
end
zeroIndex=[];
for i=1:length(step2OARs)
    if length(voxelC{step2OARs(i)})==0
        fprintf('Step 2 organ %d has zero voxel',step2OARs(i));
        return;
        
    end
end
zeroIndex=[];
for i=1:length(step3OARs)
    if length(voxelC{step3OARs(i)})==0
        planC{indexS.structures}(step3OARs(i)).structureName
        fprintf('Step 3 organ %d has zero voxel',step3OARs(i));
        return;
    end
end
choppedIM = size(planC{planC{end}.IM},2);
influenceM = getGlobalInfluenceM2(planC{indexS.IM}(choppedIM).IMDosimetry, unique([targets step2OARs step3OARs]));

numPBs = size(influenceM,2);

numPBs
sampleRateMaster = max(sampleRates);

clear struc mask3M

% Number of targets
numTargets = size(targets,2);
% create default constraints / problem setup for MOSEK
clear probDefault
probDefault.blc = [];
probDefault.buc = [];
probDefault.blx = zeros(numPBs,1);
probDefault.bux = upperBound*ones(numPBs,1);
%probDefault.ints.sub = 1:numPBs; % for integer programming


probDefault.a=[];

probDefault.qcsubi=[];
probDefault.qcsubj=[];
probDefault.qcsubk=[];
probDefault.qcval=[];

%fill a, blc, buc with hard constraints for step 1
for oar = 1:size(step1OARs,2)
    
    probDefault.a = [probDefault.a; influenceM(voxelC{step1OARs(oar)},:)];
    probDefault.blc(end+1:end+numVoxels(step1OARs(oar)),1)  = 0;
    probDefault.buc(end+1:end+numVoxels(step1OARs(oar)),1)  = DMax;
end

step = startstep-1;

fid = fopen('mosekOutput.txt','at+');
total=0;
minDoseI = [];
while step < prioptOptions.prioptCode.maxstep
    step = step+1;
    
    fprintf('\n\n\n\nStarting step %i for case %i \n',step,patID);
    
    
    % 1) prepare problem
    clear prob
    prob = probDefault; % get default constraints from previous steps
    if isfield(prob,'c')
        linOptV = prob.c;
    else
        linOptV = [];
    end
    hessianM = [];
    
    switch step
        case 1
            % setup objective function
            hessianM = sparse(numPBs,numPBs);
            linOptV  = zeros(numPBs,1);
            prob.cfix = 0; % fixed term in objective function
            
            
            for tar = 1:numTargets
                
                
                hessianM = hessianM + ( weightsPrescription(tar)/numVoxels(targets(tar)) ) * ...
                    influenceM(voxelC{targets(tar)},:)' * influenceM(voxelC{targets(tar)},:);
                linOptV = linOptV - ( weightsPrescription(tar)/numVoxels(targets(tar)) )  * ...
                    Dpre(tar) * sum(influenceM(voxelC{targets(tar)},:))';
                prob.cfix = prob.cfix + 0.5 * ( weightsPrescription(tar) * Dpre(tar) ).^2 ; % fixed term in objective function
            end
            
            if prioptOptions.prioptCode.useMaximizationOfMinDose % expand influence matrices and hessian for new variables
                
                prob.a(:,size(prob.a,2)+1 : size(prob.a,2)+numTargets) = 0;
                for tar = 1:numTargets
                    
                    aTemp = [influenceM(voxelC{targets(tar)},:) sparse(numVoxels(targets(tar)),numTargets)];
                    aTemp(:,end-numTargets+tar) = 1;
                    prob.a = [prob.a; aTemp];
                    prob.blc(end+1:end+numVoxels(targets(tar)),1) = Dpre(tar);
                    prob.buc(end+1:end+numVoxels(targets(tar)),1) = inf;
                    clear aTemp
                end
                prob.blx(end+1:end+numTargets) = 0.05 * Dpre;  % minimum dose is pushed to 95% of prescribed dose
                prob.bux(end+1:end+numTargets) = inf;
                
                hessianM = [hessianM                  sparse(numPBs,numTargets)         ; ...
                    sparse(numTargets,numPBs) speye(numTargets) ];
                linOptV  = [linOptV; zeros(numTargets,1)];
            end
            
            
        case 2
            % apply slip:
            for tar = 1:numTargets
                quadraticConstraintNew = (1+slip.quadSlipStep2) - (fixTermObjectiveFuncStep1(tar)/maxObjectiveValueStep1Abs(tar));
                prob.buc(tar) = quadraticConstraintNew;
                
            end
            
            switch lower(prioptOptions.prioptCode.csObjType)
                case 'meandose'
                    % add oars to objective function
                    linOptV  = zeros(numPBs,1);
                    for oar = 1:size(step2OARs,2)
                        
                        subObj = sum(influenceM(voxelC{step2OARs(oar)},:))' / numVoxels(step2OARs(oar)); %is mean dose when multiplied by wVector
                        linOptV = linOptV + step2Weights(oar) * subObj;
                    end
                case 'mdx'
                    %calculate length of vector x:
                    lengthV = numPBs;
                    for oar = 1:size(step2OARs,2)
                        lengthV = lengthV + numVoxels(step2OARs(oar)); %for zjs
                        lengthV = lengthV + (1+numVoxels(step2OARs(oar))) * size(alpha{oar},2); %for yas and pajs
                    end
                    
                    %calculate height of matrix A:
                    space = 0;
                    for oar = 1:size(step2OARs,2)
                        space = space + numVoxels(step2OARs(oar));
                        subspace = size(alpha{oar},2); %num alphas for this struct
                        subspace = subspace * (numVoxels(step2OARs(oar)));
                        space = space + subspace;
                    end
                    numOrigConstrs = size(prob.a,1);
                    lastAEnd = numOrigConstrs;
                    lengthA = lastAEnd+space;
                    
                    lastEnd = numPBs;
                    
                    %A matrix
                    prob.a = [prob.a sparse(numOrigConstrs,lengthV-numPBs)];
                    
                    %c vector
                    linOptV = zeros(lengthV,1); %preallocated, non-sparse
                    
                    %bx constraints
                    prob.blx = [prob.blx; sparse(lengthV-numPBs,1)];
                    prob.bux = [prob.bux; ones(lengthV-numPBs,1)*inf]; %that's all for this one
                    
                    %bc constraints
                    prob.blc = [prob.blc; sparse(lengthA-numOrigConstrs,1)]; %that's all for this one
                    prob.buc = [prob.buc; zeros(lengthA-numOrigConstrs,1)];
                    
                    disp('add oars to objective function')
                    for oar = 1:size(step2OARs,2)
                        
                        numVoxs = numVoxels(step2OARs(oar));
                        %A matrix
                        indZjsBeg = lastEnd+1;
                        indZjsEnd = lastEnd+numVoxs;
                        numVoxsOnes = ones(numVoxs,1);
                        numVoxsDiagOnes = spdiags(numVoxsOnes,[0],numVoxs,numVoxs);
                        % clear influenceM;
                        %influenceM = getSingleGlobalInfluenceM(planC{indexS.IM}(choppedIM).IMDosimetry, step2OARs(oar));
                        prob.a = [prob.a; influenceM(voxelC{step2OARs(oar)},:) ...
                            sparse(numVoxs,indZjsBeg-numPBs-1) ...
                            (-1)*numVoxsDiagOnes ...
                            sparse(numVoxs,lengthV-indZjsEnd)];
                        lastEnd = indZjsEnd;
                        
                        %c vector
                        %nothing needs to be here
                        
                        %bc constraints
                        lastAEnd = lastAEnd+numVoxs;
                        
                        %dependant on alpha:
                        for xIndex = 1:size(alpha{oar},2)
                            
                            
                            %A matrix
                            prob.a = [prob.a; sparse(numVoxs,indZjsBeg-1) ...
                                (-1)*numVoxsDiagOnes ...
                                sparse(numVoxs,lastEnd-indZjsEnd) ...
                                numVoxsOnes ...
                                numVoxsDiagOnes ...
                                sparse(numVoxs,lengthV-(lastEnd+1+numVoxs))];
                            indYas = lastEnd+1;
                            lastEnd = indYas;
                            indPajsBeg = lastEnd+1;
                            indPajsEnd = lastEnd+numVoxs;
                            lastEnd = indPajsEnd;
                            
                            %c vector
                            linOptV(indYas,1) = 1 * step2XWeights{oar}(xIndex);
                            constAVS = 1/((1-alpha{oar}(xIndex)/100)*numVoxels(step2OARs(oar)));
                            linOptV(indPajsBeg:indPajsEnd,1) = constAVS * step2XWeights{oar}(xIndex); %pajs %could also: * step2Weights(oar)
                            
                            %bx constraints
                            prob.blx(indYas,1) = -inf;
                            
                            %bc constraints
                            indBCconstrBeg = lastAEnd+1;
                            indBCconstrEnd = lastAEnd+numVoxs;
                            prob.buc(indBCconstrBeg:indBCconstrEnd,1) = inf;
                            lastAEnd = indBCconstrEnd;
                            
                        end %alpha loop
                        
                    end
            end
            disp('Done formulating Step 2.');
            
        case 3
            % apply slip:
            for tar = 1:numTargets
                quadraticConstraintNew = (1+slip.quadSlipStep3) - (fixTermObjectiveFuncStep1(tar)/maxObjectiveValueStep1Abs(tar));
                prob.buc(tar) = quadraticConstraintNew;
                
            end
            
            % reduce mean dose
            linOptV  = zeros(numPBs,1);
            for oar = 1:size(step3OARs,2)
                % clear influenceM;
                %influenceM = getSingleGlobalInfluenceM(planC{indexS.IM}(choppedIM).IMDosimetry, step3OARs(oar));
                linOptV = linOptV +  sum(influenceM(voxelC{step3OARs(oar)},:))' / numVoxels(step3OARs(oar));
            end
            linOptV = [linOptV; sparse(size(probDefault.a,2)-numPBs,1)];
            
        case 4
            % apply slip:
            for tar = 1:numTargets
                quadraticConstraintNew = (1+slip.quadSlipStep4) - (fixTermObjectiveFuncStep1(tar)/maxObjectiveValueStep1Abs(tar));
                prob.buc(tar) = quadraticConstraintNew;
            end
            i=1;
            while i<=(2*numTargets)
                probDefault.blc(minDoseI(i):minDoseI(i+1),1) = probDefault.blc(minDoseI(i):minDoseI(i+1),1) * (1-slip.minSlipStep4);
                prob.blc(minDoseI(i):minDoseI(i+1),1) = probDefault.blc(minDoseI(i):minDoseI(i+1),1);
                i = i+2;
            end
            %also allow target min dose to slip:
            %probDefault.blc(minDoseI(1):minDoseI(2),1) = probDefault.blc(minDoseI(1):minDoseI(2),1) * (1-slip.minSlipStep4);
            %prob.blc(minDoseI(1):minDoseI(2),1) = probDefault.blc(minDoseI(1):minDoseI(2),1);
            
            % find voxels with high dose (above max cutoff).
            % clear influenceM;
            
            %first:  find the influence matrix for just those voxels outside the target and margin.
            %                     doseHotspot = influenceM(voxelC{caseParams.hotspot},:) * wV;
            %                     hotVoxs = find(doseHotspot >= hotCutoff); %hotCutoff is from init_pat_xx
            %                     %we may like to 3D-expand the hotVoxs here, to allow for faster convergence.
            %                     numHotVoxs = size(hotVoxs,1);
            
            
            lengthX = size(prob.a,2); %for now
            
            %minimize sum of beam weights squared
            hessianM = sparse(lengthX, lengthX);
            B = [ones(numPBs,1)*2; sparse(lengthX-numPBs,1)];
            B = B*(1/numPBs*1/numPBs);
            hessianM = spdiags(B,0,hessianM);
            
    end
    
    
    
    prob.c = linOptV;
    eigs(hessianM)
    [prob.qosubi,prob.qosubj,prob.qoval] = find(tril(hessianM));
    clear hessianM linOptV
    
    % 2) run MOSEK (www.mosek.com)
    disp('Starting optimization...')
    param = [];
    %callback.log       = 'mosekprint';
    
    callback.loghandle = fid;
    %param.MSK_IPAR_INTPNT_NUM_THREADS = 2;
    %is useful if optimization times are really long (>> 1 min)
    param.MSK_IPAR_INTPNT_MAX_ITERATIONS = 1000;%way too small.  1000; %changed from 200, 11/9/06.
    param.MSK_DPAR_INTPNT_TOL_REL_GAP=5e-1;%1e-4; % termination criteria %uncommented 11/7/06. original commented code was 1e-4
    param.MSK_IPAR_LICENSE_WAIT = 'MSK_ON'; %added this 11/12/06
    param.MSK_IPAR_LICENSE_PAUSE_TIME = 5000; %added this 11/12/06
    param.MSK_IPAR_INFEAS_REPORT_AUTO = 'MSK_ON';
    param.MSK_DPAR_INTPNT_TOL_PFEAS=1E-4;
    param.MSK_DPAR_INTPNT_TOL_DFEAS=1E-4;
    param.MSK_DPAR_INTPNT_TOL_INFEAS=1E-4;
    %param.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_CONIC';
    param.MSK_IPAR_INTPNT_SCALING ='MSK_SCALING_NONE';
%     if maxEfficiency %for type N
%         param.MSK_IPAR_SIM_MAX_ITERATIONS = 1000;
%         param.MSK_DPAR_SIMPLEX_ABS_TOL_PIV = 1e-1 ;
%         param.MSK_IPAR_BI_MAX_ITERATIONS = 1000;
%         param.MSK_IPAR_NONCONVEX_MAX_ITERATIONS = 1000;
%         
%     end
    
    %fprintf(['Number of constraints, variables: ' num2str(size(prob.a,1)) ', ' num2str(size(prob.a,2)) '\n']);
    %fprintf(fid,['Number of constraints, variables: ' num2str(size(prob.a,1)) ', ' num2str(size(prob.a,2)) '\n']);
    
    
    %mosekopt('anapro',prob,param)
    %[r,res] = mosekopt('minimize info',prob,param,callback);
  
    [r,res] = mosekopt('minimize info',prob,param);
    fprintf('\n\n\n\nOptimization done.\n')
    fprintf('Return code: %d\n',r);
    res.sol.itr.pobjval
    total = total +res.info.MSK_DINF_OPTIMIZER_TIME;
    
    if (step > 1) %for optFormulation = 'orig'
        fprintf('Lagrangian value for slip on target 1 at step %d:  %d\n', step, res.sol.itr.suc(1));
        fprintf(fid,'Lagrangian value for slip on target 1 at step %d:  %d\n', step, res.sol.itr.suc(1));
        resC{step,3} = res.sol.itr.suc(1);
    end
    
    
    if ( isfield(res,'sol') )
        wV = res.sol.itr.xx;
      
        
    else
        wV = zeros(numPBs,1); % no solution found
    end
    
    resC{step,1} = r;
    resC{step,2} = res;
    clear r res %prob
    
    wVlong = wV; %save this for creating constraint after step 2. ('orig')
    % set wV back to real length (in case additional variable were added)
    
    wV(numPBs+1:end) = [];
    
    % set beamlets with relative weight below lowerLimit to 0
    lowerAbsLimit = prioptOptions.prioptCode.lowerLimit * max(wV);
    wV(find(wV<lowerAbsLimit))=0;
    
    % save wV into planC, into both solution and solutions for now
    planC{indexS.IM}(choppedIM).IMDosimetry.solution(step).beamletWeights = wV;
    planC{indexS.IM}(choppedIM).IMDosimetry.solutions(step).beamletWeights = wV;
    tStep=0;
    
    fprintf('\nSTEP %i:\n',step);
    fprintf(fid,'\nSTEP %i:\n',step);
    %     fprintf('Time needed (sec): %g \n\n', tStep );
    %     fprintf(fid,'Time needed (sec): %g \n', tStep );
    timesV(step) = tStep;
    
    % 3) analyze data
    
    % calc min dose, max dose, and mean dose in all targets and oars
    maxOarId = max([targets step2OARs step3OARs]);
    maxDose=zeros(maxOarId,1);
    minDose=zeros(maxOarId,1);
    meanDose=zeros(maxOarId,1);
    stdDose=zeros(maxOarId,1);
    dmaxDose=zeros(maxOarId,1);
    dminDose=zeros(maxOarId,1);
    dmeanDose=zeros(maxOarId,1);
    dstdDose=zeros(maxOarId,1);
    %dlmwrite('infmat.txt',full(influenceM),'\t');
  
    %
    for tar = [targets step2OARs step3OARs]
        % clear influenceM;
        %influenceM = getSingleGlobalInfluenceM(planC{indexS.IM}(choppedIM).IMDosimetry, tar);
        doseM = influenceM(voxelC{tar},:) * wV ;
        d = influenceM(allVoxelC{tar},:)*wV;
        if(size(doseM,1) ~= 0)
            maxDose(tar)=max(doseM);
            minDose(tar)=min(doseM);
            meanDose(tar)=mean(doseM);
            stdDose(tar)=std(doseM,1);
            
            %         if tar==11 && step==1
            %             [row,col] = find(d<minDose(tar));
            %             numVoxels(tar) = numVoxels(tar)+size(row,1);
            %             voxelC{tar} = [voxelC{tar};voxelC{tar}(row)];
            %             clear row col;
            %             [row,col] = find(d>maxDose(tar));
            %             numVoxels(tar) = numVoxels(tar)+size(row,1);
            %             voxelC{tar} = [voxelC{tar};voxelC{tar}(row)];
            %         end
            
            
        end
        
        fprintf('%22s:   min: %6.2f mean: %6.2f max: %6.2f std: %6.2f SR: %g \n',planC{indexS.structures}(tar).structureName,minDose(tar),meanDose(tar),maxDose(tar),stdDose(tar),sampleRates(tar) );
        fprintf(fid,'%22s:   min: %6.2f mean: %6.2f max: %6.2f std: %6.2f SR: %g \n',planC{indexS.structures}(tar).structureName,minDose(tar),meanDose(tar),maxDose(tar),stdDose(tar),sampleRates(tar) );
        fprintf('When All voxels used\n');
        fprintf('%22s:   min: %6.2f mean: %6.2f max: %6.2f \n',planC{indexS.structures}(tar).structureName,min(d),mean(d),max(d));

        
        if ~exist('prioptStats') || (size(prioptStats,2) < step)
            prioptStats(step) = struct('name',{{planC{indexS.structures}(tar).structureName}},'structNum',[tar],'min',[minDose(tar)],'mean',[meanDose(tar)],'max',[maxDose(tar)],'std',[stdDose(tar)],'sr',[sampleRates(tar)],'resC',[],'other',struct());
        else
            prioptStats(step).name{end+1} = planC{indexS.structures}(tar).structureName;
            prioptStats(step).structNum(end+1) = tar;
            prioptStats(step).min(end+1) = minDose(tar);
            prioptStats(step).mean(end+1) = meanDose(tar);
            prioptStats(step).max(end+1) = maxDose(tar);
            prioptStats(step).std(end+1) = stdDose(tar);
            prioptStats(step).sr(end+1) = sampleRates(tar);
        end
    end
    
    prioptStats(step).resC = resC;
    
    clear doseM tar
    fprintf('\n');
    
    % save dose to CERR
    [doseName,temp] = sprintf('Step %i',step);
    d = influenceM*wV;
    
    %[dose3DM] = inflateSampledDoseM(d,targets,sampleRateMaster);
    [dose3DM] = inflateSampledDoseM(d,unique([targets step2OARs step3OARs]),max(sampleRates));
    %[dose3DM] = getIMDose(influenceM*wV/doseScale,masterStructure,sampleRateMaster);
    IM = planC{planC{end}.IM}(choppedIM).IMDosimetry;
    if ~batchMode
        showIMDose(dose3DM,doseName);
    else
        dose2CERR(dose3DM,[],doseName,'CERR test','Test PB distribution.','UniformCT',[],'no', IM.assocScanUID); %assocScanUID added for CERR3
    end
    % save dose in planC IM structure, into both solution and solutions for now
    planC{indexS.IM}(choppedIM).IMDosimetry.solution(step).doseArray = sparse(dose3DM(:));
    planC{indexS.IM}(choppedIM).IMDosimetry.solutions(step).doseArray = sparse(dose3DM(:));
    clear doseName temp dose3DM
    fname = sprintf('%scase%d_sampling_%d.mat',saveDirectory,patID,rate);
    disp('Saving data...');
    save(fname,'planC');
    disp('Done Saving');
    
    
    % 4) add achievements to default constraints
    
    
    
    switch step
        case 1
            % add quadratic constraints
            probDefault.qcsubi = [];
            probDefault.qcsubj = [];
            probDefault.qcval  = [];
            probDefault.qcsubk = [];
            probDefault.a = [sparse(numTargets,numPBs); probDefault.a];
            probDefault.blc = [-inf*ones(numTargets,1); probDefault.blc];
            probDefault.buc = [    zeros(numTargets,1); probDefault.buc];
            maxObjectiveValueStep1 = zeros(numTargets,1);
            
            for tar=1:numTargets
                %                             hessianM = ( 1/numVoxels(targets(tar)) ) * influenceM(voxelC{targets(tar)},:)' * influenceM(voxelC{targets(tar)},:);
                %                             linOptV =  - ( 1/numVoxels(targets(tar)) ) * dosePrescribedTarget(tar) * sum(influenceM(voxelC{targets(tar)},:))';
                
                numTar = size((voxelC{targets(tar)}),1);
                
                hessianM = ( 1/numTar ) * influenceM(voxelC{targets(tar)},:)' * influenceM(voxelC{targets(tar)},:);
                linOptV =  - ( 1/numTar) * Dpre(tar) * sum(influenceM(voxelC{targets(tar)},:))';
                
                maxObjectiveValueStep1(tar)    = 0.5 * wV' * hessianM * wV + linOptV' * wV;
                fixTermObjectiveFuncStep1(tar) = 0.5 * Dpre(tar) * Dpre(tar);
                maxObjectiveValueStep1Abs(tar) = maxObjectiveValueStep1(tar) + fixTermObjectiveFuncStep1(tar);
                
                %reformulate quadratic and linear constraints for easy slip Lagrangian extraction:
                hessianM = hessianM/(maxObjectiveValueStep1Abs(tar));
                linOptV = linOptV/(maxObjectiveValueStep1Abs(tar));
                
                [qcsubi,qcsubj,qcval] = find(tril(hessianM));
                probDefault.qcsubi = [probDefault.qcsubi; qcsubi];
                probDefault.qcsubj = [probDefault.qcsubj; qcsubj];
                probDefault.qcval  = [probDefault.qcval;  qcval ];
                probDefault.qcsubk = [probDefault.qcsubk; tar * ones(size(qcsubi,1),1)];
                probDefault.a(tar,:) = sparse(linOptV');
                %                         probDefault.buc(tar) = maxObjectiveValueStep1(tar);
                qcsubi=[];
                qcsubj=[];
                qcval =[];
                qcsubk = [];
                hessianM =[];
                linOptV=[];
                tar = tar-1;
                
            end
            
            %record keeping
            prioptStats(step).other.maxObjectiveValueStep1 = maxObjectiveValueStep1;
            prioptStats(step).other.fixTermObjectiveFuncStep1 = fixTermObjectiveFuncStep1;
            prioptStats(step).other.maxObjectiveValueStep1Abs = maxObjectiveValueStep1Abs;
            
            
            clear hessianM linOptV
            clear qcsubi qcsubj qcval
            
            % add constraints for minimum doses
            if prioptOptions.prioptCode.useMaximizationOfMinDose
                for tar = 1:numTargets
                    probDefault.a=[probDefault.a; influenceM(voxelC{targets(tar)},:)];
                    mIndex = [size(probDefault.blc,1)];
                    probDefault.blc(end+1:end+numVoxels(targets(tar)),1) = minDose(targets(tar)) ;
                    minDoseI = [minDoseI,mIndex, size(probDefault.blc,1)];
                    % add constraints for maximum doses
                    
                        %probDefault.buc(end+1:end+numVoxels(targets(tar)),1) = maxDose(targets(tar)) * doseScale;
                        
                        probDefault.buc(end+1:end+numVoxels(targets(tar)),1) = maxDose(targets(tar)) ;
                   
                    
                end
            end
            
        case 2
            switch lower(prioptOptions.prioptCode.csObjType)
                case 'meandose'
                    % add constraints for mean dose
                    for oar = 1:size(step2OARs,2)
                        probDefault.a   = [probDefault.a; sum(influenceM(voxelC{step2OARs(oar)},:))/numVoxels(step2OARs(oar))];
                        probDefault.blc = [probDefault.blc;                        0];
                        probDefault.buc = [probDefault.buc; meanDose(step2OARs(oar))];
                    end
                case 'mdx'
                    probDefault.a = [prob.a; prob.c'];
                    probDefault.blc = [prob.blc; 0];
                    minObjValStep2 = prob.c' * wVlong;
                    minObjValStep2 =  minObjValStep2;
                    probDefault.buc = [prob.buc; minObjValStep2]; %no slip allotted
                    probDefault.blx = prob.blx;
                    probDefault.bux = prob.bux;
                    %probDefault.c = prob.c;
                   
                    %record keeping
                    prioptStats(step).other.minObjValStep2 = minObjValStep2;
                    %probDefault.a(1,:) = 0;
                    
            end
        case 3
            % add constraints for mean dose
            for oar = 1:size(step3OARs,2)
                probDefault.a   = [probDefault.a; sum(influenceM(voxelC{step3OARs(oar)},:))/numVoxels(step3OARs(oar)) sparse(1,size(probDefault.a,2)-numPBs)];
                probDefault.blc = [probDefault.blc;                        0];
                probDefault.buc = [probDefault.buc; meanDose(step3OARs(oar))];
            end
        case 4
            %do nothing
    end
end



clear prob wVlong


%==========================================================================
%
% END LOOP OVER STEPS
%
%==========================================================================

% final analysis:
% plot figures with beam weights, DVH's and fluence maps

s = sprintf('time%d_%d.txt',patID,rate);
fid = fopen(s,'w');
fprintf(fid,'%f',total);
fprintf('Total Time:%f\n',total);
fclose(fid);
fprintf('Slip factors used: %g, %g, %g, %g\n\n', slip.quadSlipStep2, slip.quadSlipStep3, slip.quadSlipStep4, slip.minSlipStep4);
doseNum = size(planC{planC{end}.dose},2);
l = length(targets)+length(step2OARs);
metrics=zeros(l,1);
k=1;
for tar = [targets ]
    
    metrics(k) = Dx(planC, tar, doseNum, 95);
    fprintf('D95 of %d:%f\n',tar,metrics(k));
    k=k+1;
end

for i=1:length(step2OARs)
    metrics(k) = MOHx(planC, step2OARs(i), doseNum, 5);
    fprintf('Moh5 of %d:%f\n',tar,metrics(k));
    k=k+1;
end
fprintf('\n\nFinished.\n\n');



clear temp   t   quadraticConstraintNew

clear fid fid2 fileName temp batchMode probDefault t tStep  quadraticConstraintNew





%% =======================================================================
function [prob, linOptV, constrInd, MLines] = addMeanTail(MOHType, objType, sumType, ...
    influenceM, voxelC, numVoxels, numPBs, ...
    prob, linOptV, currStructs, currStructWeights, ...
    currXVals, currXWeights, priorLengthOfVectorX, sumConstr, indivConstrs)
%adds MOH or MOC, constr or objective, to prob formulation.
%prob.c input is not used. Instead, linOptV is used.
%after calling, must say prob.c = linOptV.
%
%Assumes blx, bux, blc, buc are appropriately populated (not empty) for prev prob size.
%
%constrInd is the index of the constraint into buc. Only valid for types x01.
%MLines is the constraint/objective line for each MOH/MOC.

priorLengthOfConstrs = size(prob.a,1); %height of A %I can just get this from prob.a
constrInd = -1;

%begin.
indivConstrV = [indivConstrs{:}]; %flatten
if objType && ~sumType
    error('Sorry, an objective must be a sum (cannot be individual).');
end
if ~MOHType && ~objType %MOCType
    if sumConstr < 0 || sum(indivConstrV < 0)
        error(['MOC constraint values must be positive.  They will be converted to negative values inside the formulation.  sumConstr = ' num2str(sumConstr) ', indivConstrV = ' num2str(indivConstrV) '.']);
    end
end
%calculate length of vector x (num variables):
lengthV = priorLengthOfVectorX;
for struc = 1:size(currStructs,2)
    lengthV = lengthV + numVoxels(currStructs(struc)); %for zjs
    lengthV = lengthV + (1+numVoxels(currStructs(struc))) * size(currXVals{struc},2); %for yas and pajs
end

%calculate height of matrix A: (num constraints)
space = 0;
for struc = 1:size(currStructs,2)
    space = space + numVoxels(currStructs(struc));
    subspace = size(currXVals{struc},2); %num alphas for this struct
    subspace = subspace * (numVoxels(currStructs(struc)));
    space = space + subspace;
end
numOrigConstrs = priorLengthOfConstrs;
lastAEnd = numOrigConstrs;
lengthA = lastAEnd+space;

lastEnd = priorLengthOfVectorX;

%A matrix, first part
prob.a = [prob.a sparse(numOrigConstrs,lengthV-priorLengthOfVectorX)];

%c vector (taken care of later)
MLine = zeros(lengthV,1); %preallocated, non-sparse
MLinesHeight = 0;
for struc = 1:size(currStructs,2)
    MLinesHeight = MLinesHeight + size(currXVals{struc},2); %num alphas for this struct
end
MLines = zeros(lengthV,MLinesHeight);

%bx constraints
prob.blx = [prob.blx; sparse(lengthV-priorLengthOfVectorX,1)];
prob.bux = [prob.bux; ones(lengthV-priorLengthOfVectorX,1)*inf]; %that's all for this one

%bc constraints
prob.blc = [prob.blc; sparse(lengthA-numOrigConstrs,1)]; %that's all for this one
prob.buc = [prob.buc; zeros(lengthA-numOrigConstrs,1)];


disp('add structures to MOH constraints')
MLinesInd = 0;
for struc = 1:size(currStructs,2)
    struc
    numVoxs = numVoxels(currStructs(struc));
    
    %first, the A*w - z = 0 constraint.
    
    %A matrix
    indZjsBeg = lastEnd+1;
    indZjsEnd = lastEnd+numVoxs;
    numVoxsOnes = ones(numVoxs,1);
    numVoxsDiagOnes = spdiags(numVoxsOnes,[0],numVoxs,numVoxs);
    %influenceM = getSingleGlobalInfluenceM(planC{indexS.IM}(choppedIM).IMDosimetry, currStructs(struc));
    
    prob.a = [prob.a; influenceM(voxelC{currStructs(struc)},:) ...
        sparse(numVoxs,indZjsBeg-numPBs-1) ...
        (-1)*numVoxsDiagOnes ...
        sparse(numVoxs,lengthV-indZjsEnd)];
    lastEnd = indZjsEnd;
    
    %c vector
    %nothing needs to be here
    
    %bc constraints
    lastAEnd = lastAEnd+numVoxs;
    
    
    %now, the p-z+y >= 0 constraint:
    
    %dependant on alpha:
    for xIndex = 1:size(currXVals{struc},2)
        xIndex
        MLinesInd = MLinesInd+1;
        
        if MOHType
            zTerm = -1;
            yTerm = 1;
            pTerm = 1;
        else %MOCType
            zTerm = 1;
            yTerm = -1;
            pTerm = 1;
        end
        
        %A matrix
        prob.a = [prob.a; sparse(numVoxs,indZjsBeg-1) ...
            zTerm*numVoxsDiagOnes ... % -z
            sparse(numVoxs,lastEnd-indZjsEnd) ...
            yTerm*numVoxsOnes ... % +y
            pTerm*numVoxsDiagOnes ... % +p
            sparse(numVoxs,lengthV-(lastEnd+1+numVoxs))];
        indYas = lastEnd+1;
        lastEnd = indYas;
        indPajsBeg = lastEnd+1;
        indPajsEnd = lastEnd+numVoxs;
        lastEnd = indPajsEnd;
        
        %bc constraints
        indBCconstrBeg = lastAEnd+1;
        indBCconstrEnd = lastAEnd+numVoxs;
        prob.buc(indBCconstrBeg:indBCconstrEnd,1) = inf; % >=0 %change this to blc = -inf for <=0.
        lastAEnd = indBCconstrEnd;
        
        
        %create the MOH line, which can either be the objective function or a
        %constraint added to end of prob.a
        
        %c vector
        if MOHType
            ycTerm = 1;
            pcTerm = 1;
        else %MOCType
            ycTerm = -1;
            pcTerm = 1;
        end
        %if sumType %add one line for all
        MLine(indYas,1) = ycTerm * currXWeights{struc}(xIndex);
        constAVS = 1/((currXVals{struc}(xIndex)/100)*numVoxels(currStructs(struc)));
        MLine(indPajsBeg:indPajsEnd,1) = pcTerm * constAVS * currXWeights{struc}(xIndex) * currStructWeights(struc); %pajs
        %else %add individual lines
        %error('sumType = 0 is not yet implemented.');
        MLines(indYas,MLinesInd) = ycTerm * currXWeights{struc}(xIndex);
        MLines(indPajsBeg:indPajsEnd,MLinesInd) = pcTerm * constAVS * currXWeights{struc}(xIndex)* currStructWeights(struc); %pajs
        %end
        
        %bx constraints
        prob.blx(indYas,1) = -inf;
        
        
    end %alpha loop
    
end

%c vector:
if objType %if these MOH constraints are to be the only objective:
    linOptV = MLine;
else %MOH should be a constraint
    %let the previous objective stand:
    linOptV = [linOptV; zeros(lengthV-size(linOptV,1),1)];
    
    %turn MOH into a constraint.
    if sumType
        prob.a = [prob.a; MLine'];
        
        if MOHType
            prob.blc = [prob.blc; 0]; %could also be -inf, but we know this is >= 0 because MOH is positive.
            prob.buc = [prob.buc; sumConstr]; %no slip allotted
        else %MOCType
            prob.blc = [prob.blc; -inf];
            prob.buc = [prob.buc; -sumConstr]; %no slip allotted
        end
        constrInd = size(prob.buc,1);
    else %each individually
        prob.a = [prob.a; MLines'];
        if MOHType
            prob.blc = [prob.blc; sparse(MLinesHeight,1)];
            prob.buc = [prob.buc; indivConstrV']; %no slip allotted
        else %MOCType
            prob.blc = [prob.blc; -inf*ones(MLinesHeight,1)];
            prob.buc = [prob.buc; -indivConstrV']; %no slip allotted
        end
    end
    
    
end

