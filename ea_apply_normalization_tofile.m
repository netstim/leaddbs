function ea_apply_normalization_tofile(options,from,to,useinverse,interp,ref)
% this function applies lead-dbs normalizations to nifti files.
% currently just used to generate patient specific atlases,i.e., from MNI
% space to native space

if ~exist('interp', 'var')
    interp=4;
end





if ~exist('ref', 'var')
    ref='';
end

json = loadjson(options.subj.norm.log.method);

if ischar(interp)
    if strcmp(interp,'auto') % only works if one image supplied
        interp=detinterp(from,contains(json.method, 'ANTs'));
    end
end

if contains(json.method, 'ANTs')
    ea_ants_apply_transforms(options, from, to, useinverse, ref, '', interp);
elseif contains(json.method, 'FNIRT')
    if useinverse
        % In case inverse transformation using FSL, make sure the input has
        % exactly the same affine matrix as the template.
        for i=1:length(from)
            ea_imcalc(from{i}, [ea_space, options.primarytemplate, '.nii'], '', interp);
        end
    end
    ea_fsl_apply_normalization(options, from,to, useinverse, ref, '', interp);
elseif contains(json.method, 'SPM')
    for i=1:length(from)
        if strcmp(from{i}(end-2:end),'.gz')
            wasgz = 1;
            [~, fn] = fileparts(from{i});
            copyfile(from{i}, [tempdir, fn, '.gz']);
            gunzip([tempdir, fn, '.gz']);
            from{i} = [tempdir, fn];
        else
            wasgz = 0;
        end

        refIsGz = 0;

        % ATTENTION: when using PUSH method, transforming from native space
        % to template space should use the INVERSE transformation,
        % transforming from template space to native space should use
        % FORWARD transformation.
        usepush=1;

        if useinverse
            if isempty(ref)
                ref = options.subj.preopAnat.(options.subj.AnchorModality).coreg;
            elseif endsWith(ref, '.gz')
                refIsGz = 1;
                gunzip(ref);
                ref = strrep(ref, '.gz', '');
            end

            if usepush
                matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.subj.norm.transform.forwardBaseName, 'spm.nii']};
                matlabbatch{1}.spm.util.defs.out{1}.push.fnames = from(i);
                matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
                matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {fileparts(to{i})};
                matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {ref};
                matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
                matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0.5 0.5 0.5];
                matlabbatch{1}.spm.util.defs.out{1}.push.prefix = '';
            else
                matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.subj.norm.transform.inverseBaseName, 'spm.nii']};
                matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = from(i);
                matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fileparts(to{i})};
                matlabbatch{1}.spm.util.defs.out{1}.pull.interp = interp;
                matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
                matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0.5 0.5 0.5];
                matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
            end
            spm_jobman('run',{matlabbatch});
            clear matlabbatch
        else
            if isempty(ref)
                ref = [ea_space, options.primarytemplate, '.nii'];
            elseif endsWith(ref, '.gz')
                refIsGz = 1;
                gunzip(ref);
                ref = strrep(ref, '.gz', '');
            end

            if usepush
                matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.subj.norm.transform.inverseBaseName, 'spm.nii']};
                matlabbatch{1}.spm.util.defs.out{1}.push.fnames = from(i);
                matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
                matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {fileparts(to{i})};
                matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {ref};
                matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
                matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0.5 0.5 0.5];
                matlabbatch{1}.spm.util.defs.out{1}.push.prefix = '';
            else
                matlabbatch{1}.spm.util.defs.comp{1}.def = {[options.subj.norm.transform.forwardBaseName, 'spm.nii']};
                matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = from(i);
                matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fileparts(to{i})};
                matlabbatch{1}.spm.util.defs.out{1}.pull.interp = interp;
                matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
                matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0.5 0.5 0.5];
                matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
            end
            spm_jobman('run',{matlabbatch});
            clear matlabbatch
        end

        pth = fileparts(to{i});
        [~, fn, ext] = fileparts(from{i});

        if strcmp(to{i}(end-2:end),'.gz')
            to{i} = to{i}(1:end-3);
            gzip_output = 1;
        else
            gzip_output = 0;
        end

        movefile(fullfile(pth, ['sw', fn, ext]), to{i});

        if refIsGz
            gzip(ref);
            delete(ref);
        end

        if gzip_output
            gzip(to{i});
            delete(to{i});
        end

        if wasgz
            ea_delete([tempdir, fn, '.gz']);
            ea_delete([tempdir, fn]);
        end
    end
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



