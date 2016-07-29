function ea_ants_applytransforms(varargin)
% Wrapper for antsApplyTransforms in terms of reapplying normalizations to
% pre- and postop imaging.

ea_libs_helper;

options=varargin{1};

useinverse=0;
if nargin>1 % manual application
    fis=varargin{2};
    ofis=varargin{3};
    useinverse=varargin{4};
end

if nargin>4
    refim=varargin{5};
else
    refim=''; % use defaults
end

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    applyTransforms = [basedir, 'antsApplyTransforms.exe'];
else
    applyTransforms = [basedir, 'antsApplyTransforms.', computer('arch')];
end

directory=[options.root,options.patientname,filesep];
if nargin==1
    switch options.modality
        case 1 % MR
            fis={[directory,options.prefs.prenii_unnormalized]};
            try fis=[fis,[directory,options.prefs.tranii_unnormalized]]; end
            try fis=[fis,[directory,options.prefs.cornii_unnormalized]]; end
            try fis=[fis,[directory,options.prefs.sagnii_unnormalized]]; end

            ofis={[directory,options.prefs.gprenii]};
            try ofis=[ofis,[directory,options.prefs.gtranii]]; end
            try ofis=[ofis,[directory,options.prefs.gcornii]]; end
            try ofis=[ofis,[directory,options.prefs.gsagnii]]; end

            try lfis={[options.prefs.prenii]}; end
            try lfis=[lfis,[options.prefs.tranii]]; end
            try lfis=[lfis,[options.prefs.cornii]]; end
            try lfis=[lfis,[options.prefs.sagnii]]; end

            try
                if exist([directory,options.prefs.fa2anat],'file');
                    fis{5}=[directory,options.prefs.fa2anat];
                    ofis{5}=[directory,'gl',options.prefs.fa2anat];
                    lfis{5}=[directory,'l',options.prefs.fa2anat];
                end
            end
        case 2 % CT
            fis{1}=[directory,options.prefs.prenii_unnormalized];
            fis{2}=[directory,options.prefs.ctnii_coregistered];
            ofis{1}=[directory,options.prefs.gprenii];
            ofis{2}=[directory,options.prefs.gctnii];
            lfis{1}=[options.prefs.prenii];
            lfis{2}=[options.prefs.ctnii];
            if exist([directory,options.prefs.fa2anat],'file');
                fis{3}=[directory,options.prefs.fa2anat];
                ofis{3}=[directory,'gl',options.prefs.fa2anat];
                lfis{3}=[directory,'l',options.prefs.fa2anat];
            end
    end
end

for fi=1:length(fis)
    if ~exist(fis{fi},'file')   % skip if unnormalized file doesn't exist
        warning('%s not found. Skipping...\n',fis{fi});
        continue
    end
    if exist(ofis{fi},'file')   % skip already normalized file
        continue
    end
    % generate gl*.nii files
    [~,lprebase]=fileparts(options.prefs.prenii);
    cmd = [applyTransforms,' --verbose 1' ...
           ' --dimensionality 3 --float 1' ...
           ' -i ',ea_path_helper(fis{fi}), ...
           ' -o ',ea_path_helper(ofis{fi})];

       if useinverse
           if isempty(refim)
               refim=[directory,options.prefs.prenii_unnormalized];
           end

           if exist([directory,lprebase,'Composite.h5'],'file')
               tr=[' -r ',refim,...
                   ' -t [',ea_path_helper([directory,lprebase]),'InverseComposite.h5,0]'];
           else

               tr=[' -r ',refim,...
                   ' -t [',ea_path_helper([directory,lprebase]),'1InverseWarp.nii.gz,0]',...
                   ' -t [',ea_path_helper([directory,lprebase]),'0GenericAffine.mat,',num2str(useinverse),']'];
           end
       else
           if isempty(refim)
               refim=[options.earoot,'templates',filesep,'mni_hires.nii'];
           end
           if exist([directory,lprebase,'Composite.h5'],'file')
               tr=[' -r ',refim,...
                   ' -t [',ea_path_helper([directory,lprebase]),'Composite.h5,0]'];
           else
               tr=[' -r ',refim,...
                   ' -t [',ea_path_helper([directory,lprebase]),'1Warp.nii.gz,0]'...
                   ' -t [',ea_path_helper([directory,lprebase]),'0GenericAffine.mat,0]'];
           end
       end

       cmd=[cmd, tr];

    if ~ispc
        try
            system(['bash -c "', cmd, '"']);
        end
    else
        try
            system(cmd);
        end
    end
    
    % generate l*.nii files
    if nargin==1 % standard case
        try
            matlabbatch{1}.spm.util.imcalc.input = {[options.earoot,'templates',filesep,'bb.nii,1'];
                                                    [ofis{fi},',1']};
            matlabbatch{1}.spm.util.imcalc.output = lfis{fi};
            matlabbatch{1}.spm.util.imcalc.outdir = {directory};
            matlabbatch{1}.spm.util.imcalc.expression = 'i2';
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

            try spm_jobman('run',{matlabbatch}); end
            clear matlabbatch
        end
    end
end
