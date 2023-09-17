function ea_gencheckregpair(moving,fixed,outfn, options)
% function that uses FSL SLICER tool to create a checkreg figure.
arguments
    moving {mustBeFile} % registered moving image
    fixed {mustBeFile} % fixed image, edge will be generated based on it
    outfn {mustBeTextScalar} % output png file path
    options.preproc {mustBeNumericOrLogical} = true % preprocess the moving image
    options.lbound {mustBeNumeric} = 0.2 % lower bound for thresholding
    options.ubound {mustBeNumeric} = 99.7 % upper bound for thresholding
end

if options.preproc
    uuid = ea_generate_uuid;
    mov = ea_load_untouch_nii(moving);
    mov.img = double(mov.img);
    mov.img(isnan(mov.img)) = 0;
    SIX = sort(mov.img(mov.img(:)~=0),'descend');
    Ubound = SIX(round(length(SIX)*options.lbound/100));
    mov.img(mov.img>Ubound) = Ubound;
    Lbound = SIX(round(length(SIX)*options.ubound/100));
    mov.img(mov.img<Lbound) = Lbound;
    %mov.img(mov.img(:)~=0) = ea_contrast(mov.img(mov.img(:)~=0),2,0);

    mov.filetype = 16;
    mov.hdr.dime.scl_slope = 1;
    ea_save_untouch_nii(mov,[tempdir,'lead_temp',uuid,'.nii']);
    moving = [tempdir, 'lead_temp', uuid, '.nii'];

    threshParams = [' -e 0.05 -i ', num2str(Lbound), ' ', num2str(Ubound)];
else
    threshParams = '';
end

basedir = fileparts(outfn);
if ~exist(basedir,'dir')
    mkdir(basedir);
end

basedir = [ea_getearoot,'ext_libs',filesep,'fsl',filesep];
SLICER = ea_getExec([basedir, 'slicer'], escapePath = 1);

cmd = [SLICER, ...
    ' ', ea_path_helper(moving), ...
    ' ', ea_path_helper(fixed), ...
    threshParams, ...
    ' -a ' ea_path_helper(outfn)];

ea_runcmd(cmd, env='FSLOUTPUTTYPE=NIFTI');

if options.preproc
    ea_delete([tempdir,'lead_temp',uuid,'.nii']);
end
