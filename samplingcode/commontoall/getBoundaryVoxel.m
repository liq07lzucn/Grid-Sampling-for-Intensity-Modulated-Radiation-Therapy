for struc = unique([targets step2OARs step3OARs])
    %Find relevant voxel indices 
   
    if length(targets)>=2 && struc == targets(2)
        % subtract PTV1
        mask3M = getUniformStr(struc) & ~getUniformStr(targets(1));
       % downsample if necessary.
        if sampleRates(struc) > 1
           mask3M = getDown3Mask(mask3M, sampleRates(struc), 1) & mask3M;
        end
        mask3M = getSamplingMask(struc,mask3M,rate)& mask3M;
        voxelC{struc} = find(mask3M);
    elseif length(targets)>=3 && struc == targets(3)
        % subtract PTV1 and PTV2
        mask3M = getUniformStr(struc) & ~getUniformStr(targets(1)) & ~getUniformStr(targets(2));
        %downsample if necessary.
        if sampleRates(struc) > 1
           mask3M = getDown3Mask(mask3M, sampleRates(struc), 1) & mask3M;
        end
        mask3M = getSamplingMask(struc,mask3M)& mask3M;
        voxelC{struc} = find(mask3M);
    else    
        % all other structures
        mask3M = getUniformStr(struc);
        if sampleRates(struc) > 1
           mask3M = getDown3Mask(mask3M, sampleRates(struc), 1) & mask3M;
        end
        mask3M = getSamplingMask(struc,mask3M,rate)& mask3M;
        voxelC{struc} = find(mask3M);
        
    end    
    % Number of voxels for this structure
    numVoxels(struc) = size(voxelC{struc},1);
end
