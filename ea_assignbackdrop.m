function varargout=ea_assignbackdrop(bdstring,options,subpat,native)

if ~exist('subpat','var')
    subpat='Patient';
end
if ~exist('native','var')
    native=0; % default
try % options.native may not be defined
    if options.native
        native=1;
    else
        native=0;
    end
end
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
        haspostop=0; haspreop=0;
        try
            
            assignpatspecific(options); % use this as a probe to see if patient is defined.
            haspostop=1;
            options=ea_tempswitchoptstopre(options);
            assignpatspecific(options); % use this as a probe to see if patient is defined.
            haspreop=1;
        catch
            try
                options=ea_tempswitchoptstopre(options);
                assignpatspecific(options); % use this as a probe to see if patient is defined.
                haspreop=1;
            catch
                nopatientmode=1;
            end
        end
        if nopatientmode
            varargout{1}=ea_standardspacelist;
        else
            if native
                varargout{1}=[ea_checkhas({[subpat,' Pre-OP']},haspreop),...
                    ea_checkhas({[subpat,' Post-OP']},haspostop)];
            else
                varargout{1}=[ea_standardspacelist,...
                    ea_checkhas({[subpat,' Pre-OP']},haspreop),...
                    ea_checkhas({[subpat,' Post-OP']},haspostop)];
            end
        end


    case [subpat,' Pre-OP']
        options=ea_tempswitchoptstopre(options);
        [Vtra,Vcor,Vsag]=assignpatspecific(options);
        varargout{1}=Vtra;
        varargout{2}=Vcor;
        varargout{3}=Vsag;
    case [subpat, ' Post-OP']
        [Vtra,Vcor,Vsag]=assignpatspecific(options);
        varargout{1}=Vtra;
        varargout{2}=Vcor;
        varargout{3}=Vsag;
    case 'BigBrain 100 um ICBM 152 2009b Sym'
        if ~ea_checkinstall('bigbrain',0,1)
            ea_error('BigBrain is not installed and could not be installed automatically. Please make sure that Matlab is connected to the internet.');
        end
        varargout{1}=spm_vol(fullfile(ea_space(options),'bigbrain_2015_100um_bb.nii'));
        varargout{2}=spm_vol(fullfile(ea_space(options),'bigbrain_2015_100um_bb.nii'));
        varargout{3}=spm_vol(fullfile(ea_space(options),'bigbrain_2015_100um_bb.nii'));

    otherwise
        template=lower(strrep(bdstring,[ea_getspace,' '],''));
        varargout{1}=spm_vol(fullfile(ea_space(options),[template,'.nii']));
        varargout{2}=spm_vol(fullfile(ea_space(options),[template,'.nii']));
        varargout{3}=spm_vol(fullfile(ea_space(options),[template,'.nii']));
end


function cells=ea_checkhas(cells,has)

if ~has
    cells=cell(0);
end

function [Vtra,Vcor,Vsag]=assignpatspecific(options)
if options.native
    switch options.modality
        case 1 % MR
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii_unnormalized));
            try
                Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.cornii_unnormalized));
            catch
                Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii_unnormalized));
            end
            try
                Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.sagnii_unnormalized));
            catch
                Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii_unnormalized));
            end
        case 2 % CT
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.ctnii_coregistered));
            Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.ctnii_coregistered));
            Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.ctnii_coregistered));
            tracorpresent(1:3)=1;
    end

else
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
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tp_gctnii));
            Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tp_gctnii));
            Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tp_gctnii));
            tracorpresent(1:3)=1;
    end
end

function standardlist=ea_standardspacelist

spacedef=ea_getspacedef;
standardlist=cell(1,length(spacedef.templates));
for t=1:length(spacedef.templates)
   standardlist{t}=[spacedef.name,' ',upper(spacedef.templates{t})];
end
if strcmp(ea_getspace,'MNI_ICBM_2009b_NLIN_ASYM')
   standardlist{t+1}='BigBrain 100 um ICBM 152 2009b Sym';
end

function options=ea_tempswitchoptstopre(options)
% this generates a very temporary fake options struct that points to preop
% data instead of postop data.
if options.native
    options.prefs.tranii_unnormalized=options.prefs.prenii_unnormalized;
    options.prefs.cornii_unnormalized=options.prefs.prenii_unnormalized;
    options.prefs.sagnii_unnormalized=options.prefs.prenii_unnormalized;
    options.prefs.tp_ctnii_coregistered=options.prefs.prenii_unnormalized;
else
    options.prefs.gtranii=options.prefs.gprenii;
    options.prefs.tranii=options.prefs.prenii;
    options.prefs.gcornii=options.prefs.gprenii;
    options.prefs.cornii=options.prefs.prenii;
    options.prefs.gsagnii=options.prefs.gprenii;
    options.prefs.sagnii=options.prefs.prenii;
    options.prefs.tp_gctnii=options.prefs.gprenii;
end
