function ea_discfibers_roi_nimage_sel(obj)
% subroutine to load in rois and create an N-image for Lead-DBS
% fiberfiltering projects.

if isempty(obj.roidata)
    ea_discfibers_roi_collect(obj);
end

tdir=ea_getleadtempdir;

for side=1:size(obj.roidata.nii,2)
    cnt=1;
    for v=obj.patientselection
        obj.roidata.nii{v,side}.fname=fullfile(tdir,[ea_generate_uuid,'.nii']);
        ea_write_nii(obj.roidata.nii{v,side});
        vatlist{cnt,side}=obj.roidata.nii{v,side}.fname;
        cnt=cnt+1;
    end
end

if length(obj.patientselection)==1
    for side=1:size(vatlist,2)
        obj.roidata.nimage_sel{side}=obj.roidata.nii{obj.patientselection,side};
    end
elseif isempty(obj.patientselection)
    for side=1:size(vatlist,2)
        obj.roidata.nimage_sel{side}=struct;
    end
else
    % create N-image as well (0.5 mm isotropic)
    vatlist=[repmat({[ea_space,'t1.nii']},1,size(vatlist,2));...
        vatlist];
    for side=1:size(vatlist,2)
        matlabbatch{1}.spm.util.imcalc.input = vatlist(:,side);
        matlabbatch{1}.spm.util.imcalc.output = 'N.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {tdir};
        matlabbatch{1}.spm.util.imcalc.expression = 'ea_nansum(X(2:end,:))';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',{matlabbatch});
        clear matlabbatch
        ea_crop_nii(fullfile(tdir,'N.nii'));
        obj.roidata.nimage_sel{side}=ea_load_nii(fullfile(tdir,'N.nii'));
    end
end






