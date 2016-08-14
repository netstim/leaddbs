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

if nargin>5
    transformfile=varargin{6};
end

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    applyTransforms = [basedir, 'antsApplyTransforms.exe'];
else
    applyTransforms = [basedir, 'antsApplyTransforms.', computer('arch')];
end
directory=[options.root,options.patientname,filesep];
warpsuffix=ea_conv_antswarps(directory);

if nargin==1
    switch options.modality
        case 1 % MR
            fis={[directory,options.prefs.prenii_unnormalized]};
            ofis={[directory,options.prefs.gprenii]};
            if isfield(options.prefs,'prenii')
                lfis={[options.prefs.prenii]};
            end
            
            if isfield(options.prefs,'tranii_unnormalized')
                fis=[fis,[directory,options.prefs.tranii_unnormalized]];
                ofis=[ofis,[directory,options.prefs.gtranii]];
                lfis=[lfis,[options.prefs.tranii]];
            end
            
            if isfield(options.prefs,'cornii_unnormalized')
                fis=[fis,[directory,options.prefs.cornii_unnormalized]];
                ofis=[ofis,[directory,options.prefs.gcornii]];
                lfis=[lfis,[options.prefs.cornii]];
            end
            
            if isfield(options.prefs,'sagnii_unnormalized')
                fis=[fis,[directory,options.prefs.sagnii_unnormalized]];
                ofis=[ofis,[directory,options.prefs.gsagnii]];
                lfis=[lfis,[options.prefs.sagnii]];
            end
            if isfield(options.prefs,'fa2anat')
                if exist([directory,options.prefs.fa2anat],'file');
                    fis=[fis,[directory,options.prefs.fa2anat]];
                    ofis=[ofis,[directory,'gl',options.prefs.fa2anat]];
                    lfis=[lfis,[directory,'l',options.prefs.fa2anat]];
                end
            end
        case 2 % CT
            fis{1}=[directory,options.prefs.prenii_unnormalized];
            fis{2}=[directory,options.prefs.ctnii_coregistered];
            ofis{1}=[directory,options.prefs.gprenii];
            ofis{2}=[directory,options.prefs.gctnii];
            lfis{1}=[directory,options.prefs.prenii];
            lfis{2}=[directory,options.prefs.ctnii];
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

    % generate gl*.nii files
    [~,glprebase]=fileparts(options.prefs.gprenii);
    % use 'gl' affix for tranforms
    subdir=[options.root,options.patientname,filesep];
    try
        [~,lprebase]=fileparts(options.prefs.prenii);
%         if exist([subdir,lprebase,'Composite.h5'],'file')
%             movefile([subdir,lprebase,'Composite.h5'],[subdir,glprebase,'Composite.h5']);
%             movefile([subdir,lprebase,'InverseComposite.h5'],[subdir,glprebase,'InverseComposite.h5']);
%         end
    end
    try
        [~,lprebase]=fileparts(options.prefs.prenii);
%         if exist([subdir,lprebase,'0GenericAffine.mat'],'file')
%             movefile([subdir,lprebase,'0GenericAffine.mat'],[subdir,glprebase,'0GenericAffine.mat']);
%             try movefile([subdir,lprebase,'1Warp.nii.gz'],[subdir,glprebase,'1Warp.nii.gz']); end
%             try movefile([subdir,lprebase,'1InverseWarp.nii.gz'],[subdir,glprebase,'1InverseWarp.nii.gz']); end
%         end
    end
    
    cmd = [applyTransforms,' --verbose 1' ...
           ' --dimensionality 3 --float 1' ...
           ' -i ',ea_path_helper(fis{fi}), ...
           ' -o ',ea_path_helper(ofis{fi})];

       if useinverse
           if isempty(refim)
               refim=[directory,options.prefs.prenii_unnormalized];
           end

           if exist('transformfile','var')
               [pth,fn,ext]=fileparts(transformfile);
               tr=[' -r ',refim,...
                   ' -t [',ea_path_helper(pth),filesep,fn,ext,',0]'];
           else
                   tr=[' -r ',refim,...
                       ' -t [',ea_path_helper([directory,glprebase]),'InverseComposite',warpsuffix,',0]'];
           end
       else
           if isempty(refim)
               refim=[options.earoot,'templates',filesep,'mni_hires.nii'];
           end
           
           if exist('transformfile','var')
               [pth,fn,ext]=fileparts(transformfile);
               tr=[' -r ',refim,...
                   ' -t [',ea_path_helper(pth),filesep,fn,ext,',0]'];
           else
                   tr=[' -r ',refim,...
                       ' -t [',ea_path_helper([directory,glprebase]),'Composite',warpsuffix,',0]'];
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
end

% generate l*.nii files
if nargin==1 % standard case
    for fi=1:length(ofis)
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
