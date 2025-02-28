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

BDlist=getbdlist;

if strcmp(bdstring, 'list')
    % determine whether we are in No patient mode (could be called from
    % lead group or called from an empty patient viewer / lead anatomy
    if ~exist('options','var')
        options.patientname='';
    end

    % by default assuming patient mode, preop and postop images exist.
    nopatientmode=0;
    haspreop=1;
    haspostop=1;

    % check patient mode
    if isfield(options,'groupmode')
        nopatientmode = options.groupmode;
    elseif ~isfield(options, 'subj')
        nopatientmode=1;
    end

    % check if preop and postop images exist
    if ~nopatientmode
        if ~isfield(options.subj, 'preopAnat')
            haspreop=0;
        elseif ~isfile(options.subj.preopAnat.(options.subj.AnchorModality).coreg)
            haspreop=0;
        end

        if native
            if strcmp(options.subj.postopModality, 'MRI') && ~isfile(options.subj.postopAnat.ax_MRI.coreg)
                haspostop = 0;
            elseif strcmp(options.subj.postopModality, 'CT') && ~isfile(options.subj.postopAnat.CT.coreg)
                haspostop = 0;
            end
        else
            if strcmp(options.subj.postopModality, 'MRI') && ~isfile(options.subj.postopAnat.ax_MRI.norm)
                haspostop = 0;
            elseif strcmp(options.subj.postopModality, 'CT') && ~isfile(options.subj.postopAnat.CT.norm)
                haspostop = 0;
            end
        end
    end

    if ~haspostop && ~haspreop
        nopatientmode=1;
    end

    if nopatientmode
        varargout{1}=ea_standardspacelist;
    else
        if haspreop
            preopfiles = fieldnames(options.subj.coreg.anat.preop);
            preop=cellfun(@(x) [subpat, ' Pre-OP (' x ')'], preopfiles, 'Uniform', 0)';
        else
            preop={''};
        end

        if haspostop
            if strcmp(options.subj.postopModality, 'CT')
                postop = {[subpat,' Post-OP'], [subpat,' Post-OP (Tone-mapped)']};
            else
                postop = {[subpat,' Post-OP']};
            end
        else
            postop={''};
        end

        if native
            varargout{1}=[preop, postop];
        else
            varargout{1}=[ea_standardspacelist, preop, postop];
        end
    end

    if ~native
        % check for additional template backdrops
        for bd=1:length(BDlist{1})
            varargout{1} = [varargout{1}, BDlist{2}{bd}];
        end
    end

    % add manual choose:
    varargout{1} = [varargout{1}, {'Choose...'}];

elseif regexp(bdstring, ['^', subpat,' Pre-OP \(.*\)$'])    % pattern: "Patient Pre-OP (*)"
    whichpreop = regexp(bdstring, ['(?<=^', subpat,' Pre-OP \()(.*)(?=\))'],'match','once');
    if native
        vol = spm_vol(options.subj.preopAnat.(whichpreop).coreg);
    else
        normImage = options.subj.preopAnat.(options.subj.AnchorModality).norm;
        normImage = strrep(normImage, options.subj.AnchorModality, whichpreop);
        if ~isfile(normImage)
            ea_apply_normalization_tofile(options, options.subj.preopAnat.(whichpreop).coreg, normImage, 0, 1);
        end
    
        vol = spm_vol(normImage);
    end
    
    varargout{1} = vol;
    varargout{2} = vol;
    varargout{3} = vol;

elseif strcmp(bdstring, [subpat, ' Post-OP'])
    [Vtra,Vcor,Vsag] = assignpatspecific(options, native);

    varargout{1} = Vtra;
    varargout{2} = Vcor;
    varargout{3} = Vsag;

elseif strcmp(bdstring, [subpat, ' Post-OP (Tone-mapped)'])
    [Vtra,Vcor,Vsag] = assignpatspecific(options, native, 1);

    varargout{1} = Vtra;
    varargout{2} = Vcor;
    varargout{3} = Vsag;

elseif strcmp(bdstring, 'BigBrain 100 um ICBM 152 2009b Sym (Amunts 2013)')
%     if ~ea_checkinstall('bigbrain',0,1)
%         ea_error('BigBrain is not installed and could not be installed automatically. Please make sure that Matlab is connected to the internet.');
%     end
%     varargout{1}=spm_vol(fullfile(ea_space(options),'bigbrain_2015_100um_bb.nii'));
%     varargout{2}=varargout{1};
%     varargout{3}=varargout{1};

elseif regexp(bdstring, ['^',ea_getspace,' '])    % pattern: "MNI152NLin2009bAsym *"
    template=lower(regexp(bdstring, '(?<= )[^\W+]+(?= \()', 'match', 'once'));

    varargout{1}=spm_vol(ea_niigz(fullfile(ea_space,template)));
    varargout{2}=varargout{1};
    varargout{3}=varargout{1};
elseif strcmp(bdstring,'Choose...')
    [file,path]=uigetfile('*.nii',"MultiSelect","off");
    varargout{1}=spm_vol(ea_niigz(fullfile(path,file)));
    varargout{2}=varargout{1};
    varargout{3}=varargout{1};
elseif ismember(bdstring,BDlist{2})
    [~,ix]=ismember(bdstring,BDlist{2});
    varargout{1}=ea_load_nii([ea_space,'backdrops',filesep,BDlist{1}{ix}]);
    varargout{2}=varargout{1};
    varargout{3}=varargout{1};
else    % custom backdrop file
    varargout{1}=spm_vol(bdstring);
    varargout{2}=varargout{1};
    varargout{3}=varargout{1};
end


function [Vtra,Vcor,Vsag] = assignpatspecific(options, native, tonemapped)
scrfSuffix = '';
if isfile(options.subj.recon.recon)
    load(options.subj.recon.recon, 'reco');
    if isfield(reco, 'scrf')
        scrfSuffix = 'Scrf';
    end
end

if native
    switch options.subj.postopModality
        case 'MRI'
            if ~isempty(scrfSuffix) && ~isfile(options.subj.postopAnat.ax_MRI.coregScrf)
                ea_genscrfimages(options.subj, 'coreg');
            end
            Vtra = spm_vol(options.subj.postopAnat.ax_MRI.(['coreg', scrfSuffix]));
            if isfield(options.subj.postopAnat, 'cor_MRI') && isfile(options.subj.postopAnat.cor_MRI.(['coreg', scrfSuffix]))
                Vcor = spm_vol(options.subj.postopAnat.cor_MRI.(['coreg', scrfSuffix]));
            else
                Vcor = Vtra;
            end
            if isfield(options.subj.postopAnat, 'sag_MRI') && isfile(options.subj.postopAnat.sag_MRI.(['coreg', scrfSuffix]))
                Vsag = spm_vol(options.subj.postopAnat.sag_MRI.(['coreg', scrfSuffix]));
            else
                Vsag = Vtra;
            end
        case 'CT'
            if  ~isempty(scrfSuffix) && ~isfile(options.subj.postopAnat.CT.coregScrf)
                ea_genscrfimages(options.subj, 'coreg');
            end

            if exist('tonemapped', 'var') && tonemapped
                Vtra = spm_vol(options.subj.postopAnat.CT.(['coregTonemap', scrfSuffix]));
            else
                Vtra = spm_vol(options.subj.postopAnat.CT.(['coreg', scrfSuffix]));
            end
            Vcor = Vtra;
            Vsag = Vtra;
    end
else
    switch options.subj.postopModality
        case 'MRI'
            if  ~isempty(scrfSuffix) && ~isfile(options.subj.postopAnat.ax_MRI.normScrf)
                if ~isfile(options.subj.postopAnat.ax_MRI.coregScrf)
                    ea_genscrfimages(options.subj, 'coreg');
                end
                ea_genscrfimages(options.subj, 'norm');
            end
            Vtra = spm_vol(options.subj.postopAnat.ax_MRI.(['norm', scrfSuffix]));
            if isfield(options.subj.postopAnat, 'cor_MRI') && isfile(options.subj.postopAnat.cor_MRI.(['norm', scrfSuffix]))
                Vcor = spm_vol(options.subj.postopAnat.cor_MRI.(['norm', scrfSuffix]));
            else
                Vcor = Vtra;
            end
            if isfield(options.subj.postopAnat, 'sag_MRI') && isfile(options.subj.postopAnat.sag_MRI.(['norm', scrfSuffix]))
                Vsag = spm_vol(options.subj.postopAnat.sag_MRI.(['norm', scrfSuffix]));
            else
                Vsag = Vtra;
            end
        case 'CT'
            if  ~isempty(scrfSuffix) && ~isfile(options.subj.postopAnat.CT.normScrf)
                if ~isfile(options.subj.postopAnat.CT.coregScrf)
                    ea_genscrfimages(options.subj, 'coreg');
                end
                ea_genscrfimages(options.subj, 'norm');
            end

            if exist('tonemapped', 'var') && tonemapped
                Vtra = spm_vol(options.subj.postopAnat.CT.(['normTonemap', scrfSuffix]));
            else
                Vtra = spm_vol(options.subj.postopAnat.CT.(['norm', scrfSuffix]));
            end
            Vcor = Vtra;
            Vsag = Vtra;
    end
end


function standardlist=ea_standardspacelist
spacedef=ea_getspacedef;
standardlist=cell(1,length(spacedef.templates));
for t=1:length(spacedef.templates)
	standardlist{t}=[spacedef.name,' ',upper(spacedef.templates{t}),' (',spacedef.citation{1},')'];
end


function BDlist=getbdlist
BDlist=[{[]},{[]}]; % empty.
if exist([ea_space,'backdrops',filesep,'backdrops.txt'],'file')
    fid=fopen([ea_space,'backdrops',filesep,'backdrops.txt']);
    BDlist=textscan(fid,'%s %s');
    BDlist{2}=ea_underscore2space(BDlist{2});
end
