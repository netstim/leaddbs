function ea_apply_normalization_tofile(options,from,to,useinverse,interp,ref)
% this function applies lead-dbs normalizations to nifti files.
% currently just used to generate patient specific atlases,i.e., from MNI
% space to native space

if ischar(from)
    from = {from};
end

if ischar(to)
    to = {to};
end

if ~exist('interp', 'var')
    interp=4;
end

if ~exist('ref', 'var')
    ref='';
end

json = loadjson(options.subj.norm.log.method);

if ischar(interp)
    if strcmp(interp,'auto') % only works if one image supplied
        interp=detinterp(from,contains(json.method, {'ANTs', 'EasyReg', 'SynthMorph', 'SPM'}));
    end
end

if contains(json.method, {'ANTs', 'EasyReg', 'SynthMorph', 'SPM'})
    ea_ants_apply_transforms(options, from, to, useinverse, ref, '', interp);
elseif contains(json.method, 'FNIRT')
    if useinverse
        % In case inverse transformation using FSL, make sure the input has
        % exactly the same affine matrix as the template.
        for i=1:length(from)
            ea_imcalc({[ea_space, options.primarytemplate, '.nii'], from{i}}, from{i}, interp=interp);
        end
    end
    ea_fsl_apply_normalization(options, from,to, useinverse, ref, '', interp);
elseif contains(json.method, 'SPM')
    % Convert SPM deformation field to ITK format when necessary
    ea_convert_spm_warps(options.subj);
    ea_ants_apply_transforms(options, from, to, useinverse, ref, '', interp);
end


function interp=detinterp(from,ants)
if length(from)>1
    ea_error('Auto detection only implemented for single images.')
end
nii=ea_load_nii(from{1});
outs=unique(nii.img(:));
if ants
    interp='LanczosWindowedSinc'; % default
else
    interp=1;
end
if length(outs)<3
    if ants
        interp='GenericLabel';
    else
        interp=0;
    end
end
if length(outs)<100 % likely labeling file
    if outs==round(outs) % integers only, pretty much certainly labeling file
        if ants
            interp='GenericLabel';
        else
            interp=0;
        end
    end
end



