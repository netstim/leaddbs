function [coords_mm,trajectory,markers] = ea_runtraccore(options)

if isfile(options.subj.recon.recon)
    options.native = 0; % Load MNI reco
    [coords_mm,trajectory,markers] = ea_load_reconstruction(options);
end

if isempty(options.sides)
    return
end

% build lfile
fis={[ea_space,'bb.nii']};
switch options.subj.postopModality
    case 'MRI'
        if isfile(options.subj.norm.anat.postop.ax_MRI)
            fis = [fis; {options.subj.norm.anat.postop.ax_MRI}];
        end
        if isfield(options.subj.norm.anat.postop,'cor_MRI') && isfile(options.subj.norm.anat.postop.cor_MRI)
            fis = [fis; {options.subj.norm.anat.postop.cor_MRI}];
        end
    case 'CT'
        fis = [fis; {options.subj.norm.anat.postop.CT}];
end

tmp_dir = fullfile(options.subj.normDir,'tmp');
ea_mkdir(tmp_dir);
lpost_path = fullfile(tmp_dir, 'lpost.nii');

matlabbatch{1}.spm.util.imcalc.input = fis;
matlabbatch{1}.spm.util.imcalc.output = lpost_path;
matlabbatch{1}.spm.util.imcalc.outdir = {tmp_dir};
matlabbatch{1}.spm.util.imcalc.expression = 'sum(X(2:end,:),1)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = -1;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
spm_jobman('run',{matlabbatch});
clear matlabbatch
lnii = ea_load_untouch_nii(lpost_path);

for side = options.sides
    % call main routine reconstructing trajectory for one side.
    [coords,trajvector{side},trajectory{side},tramat] = ea_reconstruct(options,side,lnii);

    % refit electrodes starting from first electrode (this is redundant at this point).
    coords_mm{side} = ea_map_coords(coords', lpost_path)';

    [~,distmm] = ea_calc_distance(options.elspec.eldist,trajvector{side},tramat(1:3,1:3),lpost_path);

    comp = ea_map_coords([0,0,0;trajvector{side}]', lpost_path)'; % (XYZ_mm unaltered)

    trajvector{side} = diff(comp);

    normtrajvector{side} = trajvector{side}./norm(trajvector{side});

    for electrode=2:4
        coords_mm{side}(electrode,:) = coords_mm{side}(1,:)-normtrajvector{side}.*((electrode-1)*distmm);
    end

    markers(side).head = coords_mm{side}(1,:);
    markers(side).tail = coords_mm{side}(4,:);

    [xunitv, yunitv] = ea_calcxy(markers(side).head, markers(side).tail);
    markers(side).x = coords_mm{side}(1,:) + xunitv*(options.elspec.lead_diameter/2);
    markers(side).y = coords_mm{side}(1,:) + yunitv*(options.elspec.lead_diameter/2);

    coords_mm = ea_resolvecoords(markers,options);
end

% transform trajectory to mm space:
for side = options.sides
    try
        if ~isempty(trajectory{side})
            trajectory{side} = ea_map_coords(trajectory{side}', lpost_path)';
        end
    end
end

options.hybridsave = 1;
rmdir(tmp_dir, 's');
ea_methods(options,...
    ['DBS-Electrodes were automatically pre-localized in native & template space using Lead-DBS software',...
    ' (Horn & Kuehn 2015; SCR_002915; https://www.lead-dbs.org).'],...
    {'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'});
