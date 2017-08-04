function varargout=ea_assignbackdrop(bdstring,options,subpat,native)

if ~exist('subpat','var')
    subpat='Patient';
end

if ~exist('native','var')
    if isfield(options,'native')
        native = options.native;
    else
        native = 0;   % default
    end
end

whichpreop='';
if length(bdstring)>length([subpat,' Pre-OP']) && strcmp(bdstring(1:length([subpat,' Pre-OP'])),[subpat,' Pre-OP'])
    whichpreop=bdstring(length([subpat,' Pre-OP'])+1:end);
    whichpreop=upper(strrep(whichpreop,'(',''));
    whichpreop=strrep(whichpreop,')','');
    whichpreop=strrep(whichpreop,' ','');
    bdstring=[subpat,' Pre-OP'];
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
        if ~isempty(dir([options.root,options.patientname,filesep,'anat_*.nii']));
            haspreop=1;
        end
        try
            assignpatspecific(options, native); % use this as a probe to see if patient is defined.
            haspostop=1;
        end
        
        if ~haspostop && ~haspreop
            nopatientmode=1;
        end
        
        if nopatientmode
            varargout{1}=ea_standardspacelist;
        else
            if haspreop
                [options,preopfiles]=ea_assignpretra(options);
                preopfiles=cellfun(@(x) strrep(x,'anat_',''),preopfiles,'Un',0);
                preopfiles=cellfun(@(x) strrep(x,'.nii',''),preopfiles,'Un',0);
                preopfiles=cellfun(@upper,preopfiles,'Un',0);
                
                preop=cellfun(@(x) horzcat(subpat,' Pre-OP (',x,')'),preopfiles,'Un',0);
                %preop = {[subpat,' Pre-OP']};
            else
                preop={''};
            end
            postop = {[subpat,' Post-OP']};
            if native
                varargout{1}=[preop(logical(haspreop)),...
                    postop(logical(haspostop)),...
                    {'Choose...'}];
            else
                varargout{1}=[ea_standardspacelist,...
                    preop',...
                    postop(logical(haspostop)),...
                    {'Choose...'}];
            end
            
        end

    case [subpat,' Pre-OP']
        options=ea_tempswitchoptstopre(options, native, whichpreop);
        [Vtra,Vcor,Vsag]=assignpatspecific(options, native);
        varargout{1}=Vtra;
        varargout{2}=Vcor;
        varargout{3}=Vsag;
        
    case [subpat, ' Post-OP']
        [Vtra,Vcor,Vsag]=assignpatspecific(options, native);
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
        if regexp(bdstring, ['^',ea_getspace,' '])    % template has the pattern of "MNI_ICBM_2009b_NLIN_ASYM *"
            template=lower(strrep(bdstring,[ea_getspace,' '],''));
            varargout{1}=spm_vol(fullfile(ea_space(options),[template,'.nii']));
        else    % custom backdrop file
            varargout{1}=spm_vol(bdstring);
        end
        varargout{2}=varargout{1};
        varargout{3}=varargout{1};
end


function [Vtra,Vcor,Vsag]=assignpatspecific(options, native)
if native
    switch options.modality
        case 1 % MR
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii_unnormalized));
            if exist(fullfile(options.root,options.prefs.patientdir,options.prefs.cornii_unnormalized), 'file')
                Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.cornii_unnormalized));
            else
                Vcor=Vtra;
            end
            if exist(fullfile(options.root,options.prefs.patientdir,options.prefs.sagnii_unnormalized),'file')
                Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.sagnii_unnormalized));
            else
                Vsag=Vtra;
            end
        case 2 % CT
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tp_ctnii_coregistered));
            Vcor=Vtra;
            Vsag=Vtra;
    end
else
    switch options.modality
        case 1 % MR
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gtranii));
            if exist(fullfile(options.root,options.prefs.patientdir,options.prefs.gcornii),'file')
                Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gcornii));
            else
                Vcor=Vtra;
            end
            if exist(fullfile(options.root,options.prefs.patientdir,options.prefs.gsagnii),'file')
                Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gsagnii));
            else
                Vsag=Vtra;
            end
        case 2 % CT
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tp_gctnii));
            Vcor=Vtra;
            Vsag=Vtra;
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


function options=ea_tempswitchoptstopre(options, native, whichpreop,list)
% this generates a very temporary fake options struct that points to preop
% data instead of postop data.

if native
    options.prefs.tranii_unnormalized=['anat_',whichpreop,'.nii'];
    options.prefs.cornii_unnormalized=['anat_',whichpreop,'.nii'];
    options.prefs.sagnii_unnormalized=['anat_',whichpreop,'.nii'];
    options.prefs.tp_ctnii_coregistered=['anat_',whichpreop,'.nii'];
else
    [options,preniis]=ea_assignpretra(options);
    if strcmpi(preniis{1}(6:7),whichpreop)
        gfi=['glanat','.nii'];
    else
        gfi=['glanat_',whichpreop,'.nii'];
    end
    options.prefs.gtranii=gfi;
    options.prefs.gcornii=gfi;
    options.prefs.gsagnii=gfi;
    options.prefs.tp_gctnii=gfi;
    if ~exist('list','var')
        if ~exist([options.root,options.patientname,filesep,gfi],'file')
            to{1}=[options.root,options.patientname,filesep,gfi];
            from{1}=[options.root,options.patientname,filesep,['anat_',lower(whichpreop),'.nii']];
            ea_apply_normalization_tofile(options,from,to,[options.root,options.patientname,filesep],0);
        end
    end
end
