function ea_fsl_gencoordwarpfile(warpResolution,warpCoefFilename, withaffwarpfile, ref)

if ~exist(withaffwarpfile,'file')
    setenv('FSLOUTPUTTYPE','NIFTI')
    verbose=1;
%     scalingFactor = round(hdr.pixdim(2:4)'./warpResolution);
%     dataSize = hdr.dim(2:4)'.*scalingFactor;
%     data = ones(dataSize,'single');
%     hdr.pixdim(2:4) = hdr.pixdim(2:4)./scalingFactor';
%     hdr.sform44 = hdr.sform44*diag([1./scalingFactor 1]);
%     hdr.qform44 = hdr.sform44*diag([1./scalingFactor 1]);
%     
%     cbiWriteNifti(inversewarpfile, data, hdr,[],[],[],verbose);
%     clear data
    
    if ispc
        FNIRTUTILS=ea_path_helper([ea_getearoot,'ext_libs',filesep,'fsl',filesep,'fnirtfileutils','.exe']);
    else
        FNIRTUTILS=ea_path_helper([ea_getearoot,'ext_libs',filesep,'fsl',filesep,'fnirtfileutils','.', computer('arch')]);
    end
    
    command =  sprintf(' --in=%s --ref=%s --out=%s  --withaff', warpCoefFilename, ref, withaffwarpfile);
    if verbose
%        fprintf('(fslApplyWarpCoords) Computing FNIRT warp fields at a resolution of %s mm:\n',mat2str(hdr.pixdim(2:4)));
 %       disp(['  ' command])
    end;
    lcmd=[FNIRTUTILS,command];
    if ~ispc
        system(['bash -c "', lcmd, '"']);
    else
        system(lcmd);
    end
end