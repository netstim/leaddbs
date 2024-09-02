function S = ea_initializeS(varargin)

preexist = 0;

if nargin == 1
    if isstruct(varargin{1})
        options = varargin{1};
        [label, preexist] = ea_detstimname(options);
    elseif ischar(varargin{1})
        label = varargin{1};
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
            handles.stimlabel.String = label;
            if isempty(handles.stimlabel.Value)
                handles.stimlabel.Value = 1;
            end
        end
    elseif isfield(options, 'gen_newstim') && options.gen_newstim==1
        label = ea_detstimname(options);
        options.gen_newstim = 0;
        if exist('handles', 'var') && isfield(handles, 'stimlabel')
            handles.stimlabel.String = label;
            if isempty(handles.stimlabel.Value)
                handles.stimlabel.Value = 1;
            end
        end
    end
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
    load([options.subj.stimDir,filesep,ea_nt(options),S.label,filesep,'sub-', options.subj.subjId, '_desc-stimparameters.mat'], 'S');
    return
end

% right sources
for source=1:4
    for k=0:7
        eval(['S.Rs',num2str(source),'.k',num2str(k),'.perc=0;']);
        eval(['S.Rs',num2str(source),'.k',num2str(k),'.pol=0;']);
        eval(['S.Rs',num2str(source),'.k',num2str(k),'.imp=1;']);
    end
    eval(['S.Rs',num2str(source),'.amp=0;']);
    eval(['S.Rs',num2str(source),'.va=1;']);
    eval(['S.Rs',num2str(source),'.case.perc=100;']);
    eval(['S.Rs',num2str(source),'.case.pol=2;']);
    eval(['S.Rs',num2str(source),'.pulseWidth=60;']);
end

% left sources
for source=1:4
    for k=8:15
        eval(['S.Ls',num2str(source),'.k',num2str(k),'.perc=0;']);
        eval(['S.Ls',num2str(source),'.k',num2str(k),'.pol=0;']);
        eval(['S.Ls',num2str(source),'.k',num2str(k),'.imp=1;']);
    end
    eval(['S.Ls',num2str(source),'.amp=0;']);
    eval(['S.Ls',num2str(source),'.va=1;']);
    eval(['S.Ls',num2str(source),'.case.perc=100;']);
    eval(['S.Ls',num2str(source),'.case.pol=2;']);
    eval(['S.Ls',num2str(source),'.pulseWidth=60;']);
end

S.active=[1,1];
S.model='SimBio/FieldTrip (see Horn 2017)';
S.monopolarmodel=0;
S.amplitude={[0,0,0,0],[0,0,0,0]};
S=ea_activecontacts(S);
