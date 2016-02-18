function fid=ea_acpc2mni(cfg,leadfig)
% leadfig can be either a handle to the lead figure with selected patient
% directories or a cell of patient directories.
% __________________________________________________________________________________
% Copyright (C) 2016 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

% fidpoints_mm=[-0.4,1.53,-2.553       % AC
%     -0.24,-26.314,-4.393            % PC
%     -0.4,1.53,20];              % Midsag
fidpoints_mm=[0.25,1.298,-5.003       % AC
    -0.188,-24.756,-2.376            % PC
    0.25,1.298,55];              % Midsag


if iscell(leadfig)
uidir=leadfig;
else
uidir=getappdata(leadfig,'uipatdir');
end

% prompt for ACPC-coordinates:


acpc=[cfg.xmm,-cfg.ymm,cfg.zmm];
if cfg.mapmethod
    [FileName,PathName] = uiputfile('ACPC2MNI_Mapping.nii','Save Mapping...');
end

leaddir=[fileparts(which('lead')),filesep];

if ~length(uidir)
ea_error('Please choose and normalize patients first.');
end

disp('*** Converting ACPC-coordinates to MNI based on normalizations in selected patients.');
%ea_dispercent(0,'Iterating through patients');
for pt=1:length(uidir)
 %   ea_dispercent(pt/length(uidir));
    directory=[uidir{pt},filesep];
    [whichnormmethod,tempfile]=ea_whichnormmethod(directory);
%     if strcmp(whichnormmethod,'ea_normalize_ants')
%             ea_error('ANTs normalization is not supported for ACPC2MNI conversion as of now.');
%     end
    
    fidpoints_vox=ea_getfidpoints(fidpoints_mm,tempfile);
    
    [~,ptname]=fileparts(uidir{pt});
    options.prefs=ea_prefs(ptname);
    
    % warp into patient space:

    try
    [fpinsub_mm] = ea_map_coords(fidpoints_vox', tempfile, [directory,'y_ea_normparams.nii'], [directory,options.prefs.prenii_unnormalized]);
    catch
        ea_error(['Please check deformation field in ',directory,'.']);
    end
    
    fpinsub_mm=fpinsub_mm';
    
    
    fid(pt).AC=fpinsub_mm(1,:);
    fid(pt).PC=fpinsub_mm(2,:);
    fid(pt).MSP=fpinsub_mm(3,:);

    % x-dimension
    A=fid(pt).MSP-fid(pt).AC;
    B=fid(pt).PC-fid(pt).AC;
    xvec=cross(A,B); %normal to given plane
    xvec=xvec/norm(xvec);
    % y-dimension (just move from ac to pc and scale by y dimension):
    yvec=(fid(pt).PC-fid(pt).AC);
    yvec=yvec/norm(yvec);
    
    % z-dimension (just move from ac to msag plane by z dimension):
    zvec=(fid(pt).MSP-fid(pt).AC);
    zvec=zvec/norm(zvec);    
    switch cfg.acmcpc
        case 1 % relative to AC:
            warpcoord_mm=fid(pt).AC+acpc(1)*xvec+acpc(2)*yvec+acpc(3)*zvec;
        case 2 % relative to midcommissural point:
            warpcoord_mm=mean([fid(pt).AC;fid(pt).PC],1)+acpc(1)*xvec+acpc(2)*yvec+acpc(3)*zvec;
        case 3 % relative to PC:
            warpcoord_mm=fid(pt).PC+acpc(1)*xvec+acpc(2)*yvec+acpc(3)*zvec;
    end
    anat=ea_load_nii([directory,options.prefs.prenii_unnormalized]);
    warpcoord_mm=[warpcoord_mm';1];
    warpcoord_vox=anat.mat\warpcoord_mm;
    warpcoord_vox=warpcoord_vox(1:3);
    fid(pt).WarpedPointNative=warpcoord_mm(1:3)';
    
    % check it inverse normparams file has correct voxel size.
    Vinv=spm_vol([directory,'y_ea_inv_normparams.nii']);
    if ~isequal(Vinv.dim,anat.dim)
                ea_redo_inv(directory,options);
    end
        % re-warp into MNI:

        [warpinmni_mm] = ea_map_coords(warpcoord_vox, tempfile, [directory,'y_ea_inv_normparams.nii'], tempfile);
    
    warppts(pt,:)=warpinmni_mm';
    fid(pt).WarpedPointMNI=warppts(pt,:);
    
    if cfg.mapmethod==2
        anat.img(:)=0;
        anat.img(round(warpcoord_vox(1)),round(warpcoord_vox(2)),round(warpcoord_vox(3)))=1;
        anat.fname=[directory,'ACPCquerypoint.nii'];
        spm_write_vol(anat,anat.img);

        % warp into nativespace
        matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_normparams.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {[directory,'ACPCquerypoint.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {directory};
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        cfg_util('run',{matlabbatch});
        clear matlabbatch
        wfis{pt}=[directory,'wACPCquerypoint.nii'];
    end
end
%ea_dispercent(1,'end');

% create clear cut version:
if cfg.mapmethod==1
    bb=ea_load_nii([leaddir,'templates',filesep,'bb.nii']);
    bb.img(:)=0;
    warppts_vox=[warppts';ones(1,size(warppts,1))];
    warppts_vox=round(bb.mat\warppts_vox);
    
    for pnt=1:size(warppts_vox,2);
        try
            bb.img(warppts_vox(1,pnt),warppts_vox(2,pnt),warppts_vox(3,pnt))=1;
        end
    end
    
    bb.fname=[PathName,FileName];
    spm_write_vol(bb,bb.img);
elseif cfg.mapmethod==2
    % create innativespacemapped files:
    
    matlabbatch{1}.spm.util.imcalc.input = wfis;
    matlabbatch{1}.spm.util.imcalc.output = [FileName];
    matlabbatch{1}.spm.util.imcalc.outdir = {PathName};
    matlabbatch{1}.spm.util.imcalc.expression = 'sum(X)';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    cfg_util('run',{matlabbatch});
    clear matlabbatch
end

if cfg.mapmethod
% smooth clear version:

matlabbatch{1}.spm.spatial.smooth.data = {[PathName,FileName,',1']};
matlabbatch{1}.spm.spatial.smooth.fwhm = [1 1 1];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
cfg_util('run',{matlabbatch});
clear matlabbatch

[pth,fn,ext]=fileparts(bb.fname);
ea_crop_nii([pth,filesep,fn,ext]);
ea_crop_nii([pth,filesep,'s',fn,ext]);
end
assignin('base','fid',fid);




function fidpoints_vox=ea_getfidpoints(fidpoints_mm,tempfile)

V=spm_vol(tempfile);
fidpoints_vox=V(1).mat\[fidpoints_mm,ones(size(fidpoints_mm,1),1)]';
fidpoints_vox=fidpoints_vox(1:3,:)';





function o=cell2acpc(acpc)

acpc=ea_strsplit(acpc{1},' ');
if length(acpc)~=3
    acpc=ea_strsplit(acpc{1},',');
    if length(acpc)~=3
        ea_error('Please enter 3 values separated by spaces or commas.');
    end
end
for dim=1:3
o(dim,1)=str2double(acpc{dim});
end
