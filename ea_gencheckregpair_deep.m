function ea_gencheckregpair_deep(moving,fixed,outfn)
% function that uses FSLs SLICER tool to create a checkreg figure.

basedir=[ea_getearoot,'ext_libs',filesep,'fsl',filesep];
if ispc
    SLICER = ea_path_helper([basedir, 'slicer.exe']);
else
    SLICER = [basedir, 'slicer.', computer('arch')];
end

uid=ea_generate_guid;

mov=ea_load_untouch_nii([moving,'.nii']);
mov.img(isnan(mov.img))=0;
SIX=sort(mov.img(mov.img(:)~=0),'descend');
Ubound=SIX(round(length(SIX)*5/100));
mov.img(mov.img>Ubound)=Ubound;
Lbound=SIX(round(length(SIX)*95/100));
mov.img(mov.img<Lbound)=Lbound;

%mov.img(mov.img(:)~=0)=ea_contrast(mov.img(mov.img(:)~=0),2,0);
ea_save_untouch_nii(mov,[tempdir,'lead_temp',uid,'.nii'])
cmd=[SLICER,' ',ea_path_helper([tempdir,'lead_temp',uid,'.nii']),' ',ea_path_helper(fixed),' -e 0.05 -i ',num2str(Lbound),' ',num2str(Ubound),' -a ',ea_path_helper(outfn)];
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


