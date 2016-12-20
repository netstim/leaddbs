function varargout=ea_assignbackdrop(bdstring,options,subpat,native)

if ~exist('subpat','var')
    subpat='Patient';
end
if ~exist('native','var')
    native=0;
end

switch bdstring
    case 'list'
        % determine whether we are in No patient mode (could be called from
        % lead group or called from an empty patient viewer / lead anatomy
        if ~exist('options','var')
           options.patientname=''; 
        end
        if isfield(options,'groupmode')
            nopatientmode=options.groupmode;
        else
            if strcmp(options.patientname,'No Patient Selected')
                nopatientmode=1;
            else
                nopatientmode=0;
            end
        end
        try 
            assignpatspecific(options); % use this as a probe to see if patient is defined.
        catch
            nopatientmode=1;
        end
        if nopatientmode
            varargout{1}={'ICBM 152 2009b NLIN Asym T2',...
                'ICBM 152 2009b NLIN Asym T1',...
                'ICBM 152 2009b NLIN Asym PD',...
                'BigBrain 100 um ICBM 152 2009b Sym'};
        else
            if native
                varargout{1}={[subpat,' Pre-OP'],...
                    [subpat,' Post-OP']};
            else
                varargout{1}={'ICBM 152 2009b NLIN Asym T2',...
                    'ICBM 152 2009b NLIN Asym T1',...
                    'ICBM 152 2009b NLIN Asym PD',...
                    'BigBrain 100 um ICBM 152 2009b Sym',...
                    [subpat,' Pre-OP'],...
                    [subpat,' Post-OP']};
            end
        end
    case [subpat,' Pre-OP']
        options.prefs.gtranii=options.prefs.gprenii;
        options.prefs.tranii=options.prefs.prenii;
        options.prefs.gcornii=options.prefs.gprenii;
        options.prefs.cornii=options.prefs.prenii;
        options.prefs.gsagnii=options.prefs.gprenii;
        options.prefs.sagnii=options.prefs.prenii;
        [Vtra,Vcor,Vsag]=assignpatspecific(options);
        varargout{1}=Vtra;
        varargout{2}=Vcor;
        varargout{3}=Vsag;
    case [subpat, ' Post-OP'];
        [Vtra,Vcor,Vsag]=assignpatspecific(options);
        varargout{1}=Vtra;
        varargout{2}=Vcor;
        varargout{3}=Vsag;
    case 'ICBM 152 2009b NLIN Asym T2'
        varargout{1}=spm_vol(fullfile(options.earoot,'templates','mni_hires_t2.nii'));
        varargout{2}=spm_vol(fullfile(options.earoot,'templates','mni_hires_t2.nii'));
        varargout{3}=spm_vol(fullfile(options.earoot,'templates','mni_hires_t2.nii'));
    case 'ICBM 152 2009b NLIN Asym T1'
        varargout{1}=spm_vol(fullfile(options.earoot,'templates','mni_hires_t1.nii'));
        varargout{2}=spm_vol(fullfile(options.earoot,'templates','mni_hires_t1.nii'));
        varargout{3}=spm_vol(fullfile(options.earoot,'templates','mni_hires_t1.nii'));
    case 'ICBM 152 2009b NLIN Asym PD'
        varargout{1}=spm_vol(fullfile(options.earoot,'templates','mni_hires_pd.nii'));
        varargout{2}=spm_vol(fullfile(options.earoot,'templates','mni_hires_pd.nii'));
        varargout{3}=spm_vol(fullfile(options.earoot,'templates','mni_hires_pd.nii'));
    case 'BigBrain 100 um ICBM 152 2009b Sym'
        if ~ea_checkinstall('bigbrain',0,0,1)
            ea_error('BigBrain is not installed and could not be installed automatically. Please make sure that Matlab is connected to the internet.');
        end
        varargout{1}=spm_vol(fullfile(options.earoot,'templates','bigbrain_2015_100um_bb.nii'));
        varargout{2}=spm_vol(fullfile(options.earoot,'templates','bigbrain_2015_100um_bb.nii'));
        varargout{3}=spm_vol(fullfile(options.earoot,'templates','bigbrain_2015_100um_bb.nii'));
end



function [Vtra,Vcor,Vsag]=assignpatspecific(options)

switch options.modality
    case 1 % MR
        Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gtranii));
        try
            Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gcornii));
        catch
            Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gtranii));
        end
        try
            Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gsagnii));
        catch
            Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gtranii));
        end
    case 2 % CT
        Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii));
        Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii));
        Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii));
        tracorpresent(1:3)=1;
end