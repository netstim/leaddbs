function ea_gencheckregpair(moving,fixed,outfn)
% function that uses FSLs SLICER tool to create a checkreg figure.

basedir=[ea_getearoot,'ext_libs',filesep,'fsl',filesep];
if ispc
    SLICER = ea_path_helper([basedir, 'slicer.exe']);
else
    SLICER = [basedir, 'slicer.', computer('arch')];
end

cmd=[SLICER,' ',ea_path_helper(moving),' ',ea_path_helper(fixed),' -a ',ea_path_helper(outfn)];
basedir=fileparts(outfn);
if ~exist(basedir,'dir')
    mkdir(basedir);
end

setenv('FSLOUTPUTTYPE','NIFTI')
    if ~ispc
        system(['bash -c "', cmd, '"']);
    else
        system(cmd);
    end


