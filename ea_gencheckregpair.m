function ea_gencheckregpair(moving,fixed,outfn)
% function that uses FSL SLICER tool to create a checkreg figure.

basedir = [ea_getearoot,'ext_libs',filesep,'fsl',filesep];
SLICER = ea_getExec([basedir, 'slicer'], escapePath = 1);


uuid = ea_generate_uuid;

mov = ea_load_untouch_nii(moving);
mov.img = double(mov.img);
mov.img(isnan(mov.img)) = 0;
SIX = sort(mov.img(mov.img(:)~=0),'descend');
Ubound = SIX(round(length(SIX)*0.2/100));
mov.img(mov.img>Ubound) = Ubound;
Lbound = SIX(round(length(SIX)*99.7/100));
mov.img(mov.img<Lbound) = Lbound;

%mov.img(mov.img(:)~=0) = ea_contrast(mov.img(mov.img(:)~=0),2,0);
mov.filetype = 16;
mov.hdr.dime.scl_slope = 1;
ea_save_untouch_nii(mov,[tempdir,'lead_temp',uuid,'.nii'])
cmd = [SLICER,' ',ea_path_helper([tempdir,'lead_temp',uuid,'.nii']),' ',ea_path_helper(fixed),' -e 0.05 -i ',num2str(Lbound),' ',num2str(Ubound),' -a ',ea_path_helper(outfn)];
basedir = fileparts(outfn);
if ~exist(basedir,'dir')
    mkdir(basedir);
end

setenv('FSLOUTPUTTYPE','NIFTI')
ea_runcmd(cmd);

ea_delete([tempdir,'lead_temp',uuid,'.nii']);
