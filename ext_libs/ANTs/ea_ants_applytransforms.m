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
elseif isunix
    applyTransforms = [basedir, 'antsApplyTransforms.', computer];
end

subdir=[options.root,options.patientname,filesep];
if nargin==1
switch options.modality
    case 1 % MR
        fis{1}=[subdir,options.prefs.prenii_unnormalized];
        fis{2}=[subdir,options.prefs.tranii_unnormalized];
        fis{3}=[subdir,options.prefs.cornii_unnormalized];
        fis{4}=[subdir,options.prefs.sagnii_unnormalized];
        ofis{1}=[subdir,options.prefs.gprenii];
        ofis{2}=[subdir,options.prefs.gtranii];
        ofis{3}=[subdir,options.prefs.gcornii];
        ofis{4}=[subdir,options.prefs.gsagnii];
        lfis{1}=[options.prefs.prenii];
        lfis{2}=[options.prefs.tranii];
        lfis{3}=[options.prefs.cornii];
        lfis{4}=[options.prefs.sagnii];
        if exist([subdir,options.prefs.fa2anat],'file');
            fis{5}=[subdir,options.prefs.fa2anat];
            ofis{5}=[subdir,'gl',options.prefs.fa2anat];
            lfis{5}=[subdir,'l',options.prefs.fa2anat];
        end   
    case 2 % CT
        fis{1}=[subdir,options.prefs.prenii_unnormalized];
        fis{2}=[subdir,options.prefs.ctnii_coregistered];
        ofis{1}=[subdir,options.prefs.gprenii];
        ofis{2}=[subdir,options.prefs.gctnii];
        lfis{1}=[options.prefs.prenii];
        lfis{2}=[options.prefs.ctnii];
        if exist([subdir,options.prefs.fa2anat],'file');
            fis{3}=[subdir,options.prefs.fa2anat];
            ofis{3}=[subdir,'gl',options.prefs.fa2anat];
            lfis{3}=[subdir,'l',options.prefs.fa2anat];
        end
end
end

for fi=1:length(fis)
    % generate gl*.nii files
    [~,lprebase]=fileparts(options.prefs.prenii);
    cmd = [applyTransforms,' --verbose 1' ...
           ' --dimensionality 3 --float 1' ...
           ' -i ',ea_path_helper(fis{fi}), ...
           ' -o ',ea_path_helper(ofis{fi})];
           
       if useinverse
           if isempty(refim)
               refim=[subdir,options.prefs.prenii_unnormalized];
           end
           
           if exist([subdir,lprebase,'Composite.h5'],'file')
               tr=[' -r ',refim,...
                   ' -t [',ea_path_helper([subdir,lprebase]),'InverseComposite.h5,0]'];
           else
               
               tr=[' -r ',refim,...
                   ' -t [',ea_path_helper([subdir,lprebase]),'1InverseWarp.nii.gz,0]',...
                   ' -t [',ea_path_helper([subdir,lprebase]),'0GenericAffine.mat,',num2str(useinverse),']'];
           end
       else
           if isempty(refim)
               refim=[options.earoot,'templates',filesep,'mni_hires.nii'];
           end
           if exist([subdir,lprebase,'Composite.h5'],'file')
               tr=[' -r ',refim,...
                   ' -t [',ea_path_helper([subdir,lprebase]),'Composite.h5,0]'];
           else
               tr=[' -r ',refim,...
                   ' -t [',ea_path_helper([subdir,lprebase]),'1Warp.nii.gz,0]'...
                   ' -t [',ea_path_helper([subdir,lprebase]),'0GenericAffine.mat,0]'];
           end
       end
       
       cmd=[cmd,...
           tr];
       
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
        matlabbatch{1}.spm.util.imcalc.input = {[options.earoot,'templates',filesep,'bb.nii,1'];
            [ofis{fi},',1']
            };
        matlabbatch{1}.spm.util.imcalc.output = lfis{fi};
        matlabbatch{1}.spm.util.imcalc.outdir = {subdir};
        matlabbatch{1}.spm.util.imcalc.expression = 'i2';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        
         try   cfg_util('run',{matlabbatch}); end
        clear matlabbatch
    end
end
