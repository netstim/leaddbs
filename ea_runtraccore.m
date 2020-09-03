function [coords_mm,trajectory,markers]=ea_runtraccore(options)

directory = [options.root,options.patientname,filesep];

if exist([options.root,options.patientname,filesep,'ea_reconstruction.mat'],'file')
   [coords_mm,trajectory,markers]=ea_load_reconstruction(options);
end
if isempty(options.sides)
    return
end
% build lfile

fis={[ea_space,'bb.nii']};
switch options.modality
    case 1 % MR
        fis=[fis;{[directory,options.prefs.gtranii]}];
        if exist([directory,options.prefs.gcornii],'file')
            fis=[fis;{[directory,options.prefs.gcornii]}];
        end
    case 2 % CT
        fis=[fis;{[directory,options.prefs.gctnii]}];
end

% if exist([directory,options.prefs.gsagnii],'file')
%     fis=[fis;{[directory,options.prefs.gsagnii]}];
% end
matlabbatch{1}.spm.util.imcalc.input = fis;
matlabbatch{1}.spm.util.imcalc.output = [directory,'lpost.nii'];
matlabbatch{1}.spm.util.imcalc.outdir = {directory};
matlabbatch{1}.spm.util.imcalc.expression = 'sum(X(2:end,:),1)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = -1;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
spm_jobman('run',{matlabbatch});
clear matlabbatch
lnii=ea_load_untouch_nii([directory,'lpost.nii']);

for side=options.sides

    %try
    % call main routine reconstructing trajectory for one side.
    [coords,trajvector{side},trajectory{side},tramat]=ea_reconstruct(options.patientname,options,side,lnii);

    % refit electrodes starting from first electrode (this is redundant at this point).
    coords_mm{side} = ea_map_coords(coords', [directory,'lpost.nii'])';

    [~,distmm]=ea_calc_distance(options.elspec.eldist,trajvector{side},tramat(1:3,1:3),[directory,'lpost.nii']);

    comp = ea_map_coords([0,0,0;trajvector{side}]', [directory,'lpost.nii'])'; % (XYZ_mm unaltered)

    trajvector{side}=diff(comp);

    normtrajvector{side}=trajvector{side}./norm(trajvector{side});

    for electrode=2:4
        coords_mm{side}(electrode,:)=coords_mm{side}(1,:)-normtrajvector{side}.*((electrode-1)*distmm);
    end

    markers(side).head=coords_mm{side}(1,:);
    markers(side).tail=coords_mm{side}(4,:);

    orth=null(normtrajvector{side})*(options.elspec.lead_diameter/2);

    markers(side).x=coords_mm{side}(1,:)+orth(:,1)';
    markers(side).y=coords_mm{side}(1,:)+orth(:,2)'; % corresponding points in reality

    coords_mm=ea_resolvecoords(markers,options);
end

% transform trajectory to mm space:
for side=1:length(options.sides)
    try
        if ~isempty(trajectory{side})
            trajectory{side}=ea_map_coords(trajectory{side}', [directory,'lpost.nii'])';
        end

    end
end

options.hybridsave=1;
ea_delete([directory,'lpost.nii']);
ea_methods(options,...
    ['DBS-Electrodes were automatically pre-localized in native & template space using Lead-DBS software',...
    ' (Horn & Kuehn 2015; SCR_002915; https://www.lead-dbs.org).'],...
    {'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'});
