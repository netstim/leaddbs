function [AllX,space]=ea_exportefieldmapping(vatlist,obj)

disp('Need to export Efields in proper format, this may take a while');
outdir=[fileparts(obj.leadgroup),filesep,'sweetspots',filesep,obj.ID,filesep];
ea_mkdir(outdir);
copyfile([ea_space,'bb.nii'],[outdir,'bb_nan.nii']);
nii=ea_load_nii([outdir,'bb_nan.nii']);
nii.dt(1) = 16;
nii.img(:)=nan;
ea_write_nii(nii);
allV{1}=[outdir,'bb_nan.nii'];

cnt=2;
for vat=1:size(vatlist,1)
    for side=1:size(vatlist,2)
        allV{cnt}=vatlist{vat,side};
        cnt=cnt+1;
    end
end

% export mean to get bounding box
matlabbatch{1}.spm.util.imcalc.input = allV';
matlabbatch{1}.spm.util.imcalc.output = 'efield_bb.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {outdir};
matlabbatch{1}.spm.util.imcalc.expression = 'ea_nansum(X)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 16;

spm_jobman('run',{matlabbatch});
clear matlabbatch

if size(vatlist,2)>1
    sidesuffices={'_r','_l'};
    ea_split_nii_lr([outdir,'efield_bb.nii']);
else
    sidesuffices={''};
end

% create final space:
for side=1:size(vatlist,2)
    nii=ea_load_nii([outdir,'efield_bb',sidesuffices{side},'.nii']);
    if ~isfield(obj.M,'pseudoM')
        nii.img(nii.img<150)=nan;
    end
    ea_write_nii(nii);
    ea_crop_nii(nii.fname);
    nii=ea_load_nii(nii.fname); % reload for space function.
    space{side}=nii;
end

% now conform each VTA to space
AllX=cell(size(vatlist,2),1);
for vat=1:size(vatlist,1)
    for side=1:size(vatlist,2)
        copyfile(vatlist{vat,side},[outdir,'tmp_efield.nii']);
        ea_conformspaceto([outdir,'efield_bb',sidesuffices{side},'.nii'],...
            [outdir,'tmp_efield.nii'],0);
        nii=ea_load_nii([outdir,'tmp_efield.nii']);

       if ~exist('AllX','var')
          AllX{side}=zeros(size(vatlist,1),numel(nii.img));
       end
       AllX{side}(vat,:)=nii.img(:);
    end
end
ea_delete([outdir,'efield_bb.nii']);
ea_delete([outdir,'efield_bb_l.nii']);
ea_delete([outdir,'efield_bb_r.nii']);
ea_delete([outdir,'tmp_efield.nii']);
ea_delete([outdir,'bb_nan.nii']);
ea_delete(outdir);
