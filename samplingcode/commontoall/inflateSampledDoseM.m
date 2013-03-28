function [ dose3D ] = inflateSampledDoseM( doseInStruct,structsV,sampleRate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    global planC
    indexS = planC{end};
    dose3D = [];
    scanNumV = getStructureAssociatedScan(structsV);
    scanNum = scanNumV(1);
    %For each struct...
    for structNum = structsV
   
        if isempty(dose3D)
            %sV = getUniformizedSize(planC);
            sV = getUniformScanSize(planC{indexS.scan}(scanNum));
            %sV(3) = 45;
            dose3D = zeros(sV);
        end

        %If sub-sampled, use 3-D interpolation to fill out dose.
        if sampleRate ~= 1

            disp('Inflating downsampled dose distribution...')
            if rem(log2(sampleRate),1) ~= 0
                error('Sample factor must (currently) be a power of 2.')
            end

            %Get interpolation coords.
            maskSample3D = getDown3Mask(dose3D, sampleRate, 1);
            maskStruct3D = getUniformStr(structNum);
            %maskStruct3D(3) = 45;
          %  maskSample3D(3) = 45;
            [rInterpV, cInterpV, sInterpV] = find3d(maskStruct3D & ~maskSample3D);

            %Now do 3-D interpolation:
            sizeV = size(dose3D);        

            doseInterp = doseInStruct;
            whereV = find(maskSample3D(:));
            clear maskSample3D;
            doseInterp = doseInterp(whereV);
            doseInterp = reshape(doseInterp,[sizeV(1)/2^log2(sampleRate),sizeV(2)/2^log2(sampleRate),sizeV(3)]);

            rDownV = (rInterpV + 2^log2(sampleRate) - 1)/2^log2(sampleRate);
            cDownV = (cInterpV + 2^log2(sampleRate) - 1)/2^log2(sampleRate);

            fillsV = matInterp3(rDownV,cDownV,sInterpV,doseInterp);
            clear doseInterp;

            ind3V = sub2ind(sizeV,rInterpV,cInterpV,sInterpV);  %index for interpolated values
            doseInStruct(ind3V) = fillsV;
            %doseInStruct(~maskStruct3D) = 0;
            clear maskStruct3D;
            disp('Finished inflating.')
        end


        for i = 1:size((dose3D),3),
            tmp = [dose3D(:,:,i)==0];
            indStart = sub2ind(sV,1,1,i);
            indStop  = sub2ind(sV,sV(1),sV(2),i);
            doseSlice = reshape(doseInStruct(indStart:indStop),sV(1),sV(2));
            tmp1 = doseSlice .* tmp;  %voxel over-adding avoided
            tmp2 = dose3D(:,:,i) + tmp1;
            dose3D(:,:,i) = tmp2;
        end
   
    end
    clear doseInStruct;
end

