
function varargout = getUniformStrWithSampling(structNumV,rate, planC, optS, generateData)
%"getUniformStr"
%   Get a 3D representation of a structure or union of structures,
%   registered to the uniformized scan.
%
%   In the case of multiple scanSets, the structure is registered to the
%   uniformized scan of its associated scanSet.  If multiple structures are
%   requested, they must all have the same associated scanSet.
%
%   If the number of output arguments is 1 or 2, the first argument is a 3D
%   mask representing the structure(s).  See first example.
%
%   If the number of output arguments is 3 or 4, the first 3 outputs are
%   vectors listing the row, column and slice coordinates of voxels in the
%   structure(s).  See second example.
%
%   structNumV can be a vector of multiple structures, in which case the
%   union of the masks of those structures is returned.
%
%   The optional parameter generateData is set to 1 if the uniformized data
%   should be generated for structures that appear un-uniformized.  0 if
%   empty data should be returned for such structures.  1 is the default.
%
%   LM:  9 may 03, JOD, create uniformized structures if they don't exist.
%                       Convert input string to a number in case structure
%                       name was specified.
%       12 Dec 04, JRA, Now supports multiple uniform scanSets.
%
%Usage:
%   function [S, planC] = getUniformStr(structNumV, planC, optS, generateData);
%   function [r,c,s, planC] = getUniformStr(structNumV, planC, optS, generateData);
%
% Copyright 2010, Joseph O. Deasy, on behalf of the CERR development team.
%
% This file is part of The Computational Environment for Radiotherapy Research (CERR).
%
% CERR development has been led by:  Aditya Apte, Divya Khullar, James Alaly, and Joseph O. Deasy.
%
% CERR has been financially supported by the US National Institutes of Health under multiple grants.
%
% CERR is distributed under the terms of the Lesser GNU Public License.
%
%     This version of CERR is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
% CERR is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CERR.  If not, see <http://www.gnu.org/licenses/>.


%Check if plan passed, if not use global.

if ~exist('planC')
    global planC;
end
indexS = planC{end};

if ~exist('generateData')
    generateData = 1;
end

%Determine whether to return mask or r,c,s coordinates.
if nargout == 0 | nargout == 1 | nargout == 2
    wantMask = 1; %Return 3D mask
elseif nargout == 3 | nargout == 4
    wantMask = 0; %Return RCS.
else
    error('Invalid number of output arguments in call to getUniformStr.');
end

%Check if optS passed.  If not, try using global stateS options, else use
%options saved in the CERR plan.
if ~exist('optS')
    try
        global stateS;
        optS = stateS.optS;
    catch
        optS = planC{indexS.CERROptions};
    end
end

%If structNumV is a string or cell of strings, convert to numerics.
if isstr(structNumV) | iscell(structNumV)
    if isstr(structNumV)
        structNumV = {structNumV};
    end
    
    numericStructList = [];
    for i=1:length(structNumV)
        Str_Name = structNumV{i};
        structures = lower({planC{indexS.structures}.structureName});
        [jnk, ia, ib] = intersect(lower(structNumV{i}), structures);
        structNum = ib;
        if isempty(structNum)
            error(['Structure ''' Str_Name ''' does not exist in plan.']);
        end
        numericStructList = [numericStructList structNum];
    end
    structNumV = numericStructList;
end

%Determine which scanSet the structure(s) are registered to.
[scanNum, relStructNum] = getStructureAssociatedScan(structNumV, planC);
if length(unique(scanNum)) > 1
    error('All structures passed to getUniformStr must be registered to the same scan.');
end
scanNum = scanNum(1);

%Get uniformized structure data, generating if necessary.
[indicesC, bitsC, planC] = getUniformizedData(planC, scanNum, 'yes');

%Get size of uniform scan.
siz = getUniformScanSize(planC{indexS.scan}(scanNum));
numRows     = siz(1);
numCols     = siz(2);
numSlices   = siz(3);

[x,y,z] = getUniformScanXYZVals(planC{indexS.scan}(scanNum));

%Convert indicesM and bitsM into a single full array:
%use boolean to save memory.
if wantMask
    S = logical(repmat(uint8(0), siz));
end

%maskV = repmat(logical(0), [length(bitsM) 1]);
otherIndicesM = [];
for i=1:length(structNumV)
    
    structNum = structNumV(i);
    
    %Get the cell index
    if relStructNum(i)<=52
        cellNum = 1;
    else
        cellNum = ceil((structNum-52)/8)+1;
    end
    
    indicesM = indicesC{cellNum};
    bitsM    = bitsC{cellNum};
    
    %eliminate uniformized data that, for whatever reason, lies outside the
    %defined region.
    bitsM(find(indicesM(:,3)>numSlices),:) = [];
    indicesM(find(indicesM(:,3)>numSlices),:) = [];
    
    %get values corresponding to requested structure.
    if relStructNum(i)<=52
        bitMaskV = logical(bitget(bitsM(:),relStructNum(i)));
    else
        bitMaskV = logical(bitget(bitsM(:),relStructNum(i)-52-8*(cellNum-2)));
    end
    
    %If no voxels are part of structNum, suggests that the structure needs
    %to be uniformized. Uniformize it.  This is recursive so if a structure
    %exists that cannot be uniformized...trouble!
    if ~any(bitMaskV) & generateData
        warning(['Structure ' planC{indexS.structures}(structNum).structureName ' does not appear to be uniformized.  Adding it to uniformized data.']);
        planC = updateStructureMatrices(planC, structNum);
        otherIndices = [];
        [otherIndices(:,1), otherIndices(:,2), otherIndices(:,3), planC] = getUniformStr(structNumV, planC, optS, 0);
        break;
    end
    
    otherIndices = reshape(indicesM([bitMaskV bitMaskV bitMaskV]), [], 3);
    
    clear bitMaskV;
    
    indPresentV = ismember(otherIndices,otherIndicesM,'rows');
    otherIndicesM = [otherIndicesM; otherIndices(~indPresentV,:)];
    
    clear otherIndices
    
end
clear bitsM maskV indicesM;

if wantMask
    
    [boundaryIndices,rows] = getBoundaryIndices(otherIndicesM,x,y,z,structNumV);
    %otherIndicesM(rows,:)=[];
    remainingIndices = setdiff(otherIndicesM,boundaryIndices,'rows');
    %innerIndicies = getpVoxels(otherIndicesM,x,y,z,structNumV,10);
    innerIndicies = getpVoxelsGrid(remainingIndices,x,y,z,structNumV,rate);
    %innerIndicies = getpVoxelsRandom(otherIndicesM,x,y,z,structNumV,10);
    fprintf('Structure:%d Boundary:%d',structNumV,size(boundaryIndices,1));
    combinedIndices = [boundaryIndices;innerIndicies];
    fprintf('Inner:%d Total:%d\n',size(innerIndicies,1),size(otherIndicesM,1));

    
    %	boundaryIndices = otherIndicesM;
    
    %     if structNumV == 4
    %         boundaryIndices = otherIndicesM;
    %     end
    %Convert 3D indices to vectorized indices.
    if ~isempty(boundaryIndices)
        %indexV = sub2ind([numRows,numCols,numSlices],double(otherIndicesM(:,1)),double(otherIndicesM(:,2)),double(otherIndicesM(:,3)));
        indexV = sub2ind([numRows,numCols,numSlices],double(combinedIndices(:,1)),double(combinedIndices(:,2)),double(combinedIndices(:,3)));
        
    else
        indexV = [];
    end
    %    k=1;
    %     for i=1:size(indexV,1)
    %         x = indexV(i);
    %         if ismember(x,indexV1) ~= 0
    %             tmpIndex(k)=x;
    %             k=k+1;
    %         end
    %     end
    %     indexV = tmpIndex';
    %     s = size(indexV,1);
    %     indexV(1:1:3*s/4,:) = [];
    %     clear otherIndicesM;
    
    %Fill S array in with values where the structure exists.
    S(indexV) = 1;
    varargout{1} = S;
    varargout{2} = planC;
else
    %Assign r,c,s coordinates to output.
    %     varargout{1} = otherIndicesM(:,1);
    %     varargout{2} = otherIndicesM(:,2);
    %     varargout{3} = otherIndicesM(:,3);
    
    varargout{1} = boundaryIndices(:,1);
    varargout{2} = boundaryIndices(:,2);
    varargout{3} = boundaryIndices(:,3);
    clear otherIndicesM;
    varargout{4} = planC;
end
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
            
            %also allow target min dose to slip:
            probDefault.blc(minDoseI(1):minDoseI(2),1) = probDefault.blc(minDoseI(1):minDoseI(2),1) * (1-slip.minSlipStep4);
            prob.blc(minDoseI(1):minDoseI(2),1) = probDefault.blc(minDoseI(1):minDoseI(2),1);
            
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
    param.MSK_DPAR_INTPNT_TOL_REL_GAP=1e-1;%1e-4; % termination criteria %uncommented 11/7/06. original commented code was 1e-4
    param.MSK_IPAR_LICENSE_WAIT = 'MSK_ON'; %added this 11/12/06
    param.MSK_IPAR_LICENSE_PAUSE_TIME = 5000; %added this 11/12/06
    param.MSK_IPAR_INFEAS_REPORT_AUTO = 'MSK_ON';
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
                    minDoseI = [size(probDefault.blc,1)];
                    probDefault.blc(end+1:end+numVoxels(targets(tar)),1) = minDose(targets(tar)) ;
                    minDoseI = [minDoseI, size(probDefault.blc,1)];
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

s = sprintf('time%d.txt',patID);
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

           for oar = 1:size(step3OARs,2)
                            probDefault.a   = [probDefault.a; sum(influenceM(voxelC{step3OARs(oar)},:))/numVoxels(step3OARs(oar)) sparse(1,size(probDefault.a,2)-numPBs)];
                            probDefault.blc = [probDefault.blc;                        0];
                            probDefault.buc = [probDefault.buc; meanDose(step3OARs(oar))*doseScale];
                        end
                    case 4
                        %do nothing
                end
                
                
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
    if ~batchMode
        plotBeamWeights(prioptOptions.prioptCode.maxstep,beamWeightFigure);
        plotFluenceMaps(prioptOptions.prioptCode.maxstep,fluenceMapFigure);
        plotDVHs(prioptOptions.prioptCode.maxstep,prioptOptions.prioptCode.useMultiplePlots,dvhStructures,sampleRates,dvhFigure,voxelC);
    end
    tTotal = toc;
    fprintf('Rate:%d\n',rte);
    calcmohdx
    
    fprintf('Total time needed for optimization (sec): %g    (', tTotal );
    fprintf(' %5.0f',timesV);fprintf(')\n');
    fprintf(fid,'\n\nTotal time needed for optimization (sec): %g    (', tTotal );
    fprintf(fid,' %5.0f',timesV);fprintf(fid,')\n');
    fprintf('Slip factors used: %g, %g, %g, %g\n\n', slip.quadSlipStep2, slip.quadSlipStep3, slip.quadSlipStep4, slip.minSlipStep4);
    fprintf(fid,'Slip factors used: %g, %g, %g, %g\n\n', slip.quadSlipStep2, slip.quadSlipStep3, slip.quadSlipStep4, slip.minSlipStep4);
    
    
    fprintf('\n\nFinished.\n\n');
    
    
    clear temp   t quadraticConstraintNew
    diffrtTim(end+1) = tTotal/60;
    if(rte == 10)
        tm = tTotal/60;
      
    end
    %clear planC;
    disp('Loading data...');
    extension='';
    pfname = sprintf('%scase%d.mat',saveDirectory,patID);
    load(pfname);
    disp('Done Loading');
end
clear fid  fileName temp batchMode probDefault t tStep tTotal quadraticConstraintNew






%% =======================================================================




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

