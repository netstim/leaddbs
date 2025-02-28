function ea_discfibers_roi_collect(obj)
% subroutine to load in rois and create an N-image for Lead-DBS
% fiberfiltering projects.

if isempty(obj.roidata)
    if isfield(obj.M,'pseudoM')
        vatlist = obj.M.ROI.list;
    else
        vatlist = ea_discfibers_getvats(obj);
    end

    tdir=ea_getleadtempdir;

    for side=1:size(vatlist,2)
        for v=1:size(vatlist,1)
            if ~isfile(vatlist{v,side})
                obj.roidata.nii{v,side} = [];
                vatlist{v,side} = [];
                ea_cprintf('CmdWinWarnings', 'vatlist(%d,%d) not found when loading ROIs!\n', v, side);
            else
                obj.roidata.nii{v,side}=ea_load_nii(vatlist{v,side});
                if endsWith(vatlist{v,side}, '.gz')
                    obj.roidata.nii{v,side}.fname=fullfile(tdir,[ea_generate_uuid,'.nii']);
                    ea_write_nii(obj.roidata.nii{v,side});
                    vatlist{v,side}=obj.roidata.nii{v,side}.fname;
                end
            end
        end
    end

    % create N-image as well (0.5 mm isotropic)
    vatlist=[repmat({[ea_space,'t1.nii']},1,size(vatlist,2));...
        vatlist];
    for side=1:size(vatlist,2)
        fileList = vatlist(:,side);
        fileList(cellfun(@isempty, fileList)) = [];
        matlabbatch{1}.spm.util.imcalc.input = fileList;
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
        obj.roidata.nimage{side}=ea_load_nii(fullfile(tdir,'N.nii'));
    end
end






