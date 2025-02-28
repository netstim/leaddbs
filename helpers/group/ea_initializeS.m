function S = ea_initializeS(varargin)

preexist = 0;

if nargin == 1
    if ischar(varargin{1}) % only label provided
        label = varargin{1};
        options.elspec.numContacts = 4;
        ea_cprintf('CmdWinWarnings', 'Missing number of contacts (options.elspec.numContacts) when initializing "S"! Set to 4 by default.\n');
    elseif isstruct(varargin{1}) % 'options' provided
        options = varargin{1};
        [label, preexist] = ea_detstimname(options);
    end
elseif nargin > 1
    label = varargin{1};
    options = varargin{2};

    if nargin == 3
        handles = varargin{3};
    end

    if isempty(label)
        [label, preexist] = ea_detstimname(options);
        if exist('handles', 'var') && isfield(handles, 'stimlabel')
            if ~isempty(label)
                handles.stimlabel.String = label;
                if isempty(handles.stimlabel.Value)
                    handles.stimlabel.Value = 1;
                end
            else
                S = [];
                return;
            end
        end
    elseif isfield(options, 'gen_newstim') && options.gen_newstim==1
        label = ea_detstimname(options);
        options.gen_newstim = 0;
        if exist('handles', 'var')
            figure(handles.stimfig);
            if isfield(handles, 'stimlabel')
                if ~isempty(label)
                    handles.stimlabel.String = label;
                    if isempty(handles.stimlabel.Value)
                        handles.stimlabel.Value = 1;
                    end
                else
                    S = [];
                    return;
                end
            end
        end
    end
end

% options.elspec could be missing when switching simulation model in LeadGroup
if ~isfield(options, 'elspec') && strcmp(options.leadprod, 'group')
    elstruct = getappdata(handles.stimfig, 'elstruct');
    actpt = getappdata(handles.stimfig, 'actpt');
    options.elmodel = elstruct(actpt).elmodel;
    options = ea_resolve_elspec(options);
end

if isfield(options, 'UsePreExistingStim')
    preexist = options.UsePreExistingStim;
end

if ~iscell(label)
    label = {label};
end

if exist('handles', 'var') && isfield(handles, 'stimlabel') && ~isempty(handles.stimlabel.Value)
    S.label = label{handles.stimlabel.Value};
else
    S.label = label{1};
end

if preexist
    stimParamFile = [options.subj.stimDir,filesep,ea_nt(options),S.label,filesep,'sub-', options.subj.subjId, '_desc-stimparameters.mat'];
    if isfile(stimParamFile)
        S = ea_checkStimParams(stimParamFile);
        existingNumContacts = sum(startsWith(fieldnames(S.Rs1), 'k'));
        if existingNumContacts < options.elspec.numContacts
            for source=1:4
                for k=existingNumContacts+1:options.elspec.numContacts
                    S.(['Rs',num2str(source)]).(['k',num2str(k)]).perc=0;
                    S.(['Rs',num2str(source)]).(['k',num2str(k)]).pol=0;
                    S.(['Rs',num2str(source)]).(['k',num2str(k)]).imp=1;
                    S.(['Ls',num2str(source)]).(['k',num2str(k)]).perc=0;
                    S.(['Ls',num2str(source)]).(['k',num2str(k)]).pol=0;
                    S.(['Ls',num2str(source)]).(['k',num2str(k)]).imp=1;
                end
            end
            S.numContacts = options.elspec.numContacts;
        elseif existingNumContacts > options.elspec.numContacts
            for source=1:4
                excessiveFields = strcat('k', arrayfun(@num2str, options.elspec.numContacts+1:existingNumContacts, 'Uni', 0));
                S.(['Rs',num2str(source)]) = rmfield(S.(['Rs',num2str(source)]), excessiveFields);
                S.(['Ls',num2str(source)]) = rmfield(S.(['Ls',num2str(source)]), excessiveFields);
            end
            S.numContacts = options.elspec.numContacts;
        end
        return
    end
end

% right sources
for source=1:4
    for k=1:options.elspec.numContacts
        S.(['Rs',num2str(source)]).(['k',num2str(k)]).perc=0;
        S.(['Rs',num2str(source)]).(['k',num2str(k)]).pol=0;
        S.(['Rs',num2str(source)]).(['k',num2str(k)]).imp=1;
    end
    S.(['Rs',num2str(source)]).case.perc = 100;
    S.(['Rs',num2str(source)]).case.pol = 2;
    S.(['Rs',num2str(source)]).amp = 0;
    S.(['Rs',num2str(source)]).va = 2;
    S.(['Rs',num2str(source)]).pulseWidth = 60;
end

% left sources
for source=1:4
    for k=1:options.elspec.numContacts
        S.(['Ls',num2str(source)]).(['k',num2str(k)]).perc=0;
        S.(['Ls',num2str(source)]).(['k',num2str(k)]).pol=0;
        S.(['Ls',num2str(source)]).(['k',num2str(k)]).imp=1;
    end
    S.(['Ls',num2str(source)]).case.perc = 100;
    S.(['Ls',num2str(source)]).case.pol = 2;
    S.(['Ls',num2str(source)]).amp = 0;
    S.(['Ls',num2str(source)]).va = 2;
    S.(['Ls',num2str(source)]).pulseWidth = 60;
end

S.active = [1,1];
S.model = 'SimBio/FieldTrip (see Horn 2017)';
S.monopolarmodel = 0;
S.amplitude = {[0,0,0,0],[0,0,0,0]};
S.numContacts = options.elspec.numContacts;
S = ea_activecontacts(S);
S.sources = 1:4;
S.volume = [0 0];
S.ver = '2.0'; % refined version of S
