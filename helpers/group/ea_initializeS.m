function S = ea_initializeS(varargin)

if nargin>1
    options=varargin{2};
    handles=varargin{3};
end

if nargin
    if isempty(varargin{1})
        [labels, preexist] = ea_detstimname(options);
        set(handles.stimlabel, 'String', labels);
    elseif (isfield(options, 'gen_newstim') && options.gen_newstim==1)
        labels = ea_detstimname(options);
        preexist = 0;
        options.gen_newstim = 0;
        set(handles.stimlabel, 'String', labels);
    else
        labels = varargin{1};
        preexist = 0;
    end
else
    [labels, preexist] = ea_detstimname(options);
end

if ~iscell(labels)
    labels={labels};
end

try
    S.label = labels{get(handles.stimlabel,'Value')};
catch
    S.label = labels{1};
end

if preexist
    load([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),S.label,filesep,'stimparameters.mat']);
    return
end

% right sources
for source=1:4
    for k=0:7
        eval(['S.Rs',num2str(source),'.k',num2str(k),'.perc=0;']);
        eval(['S.Rs',num2str(source),'.k',num2str(k),'.pol=0;']);
        eval(['S.Rs',num2str(source),'.k',num2str(k),'.imp=1;']);
    end
    
    if contains(options.elmodel, 'Aleva')
        for k=8:11
            eval(['S.Rs',num2str(source),'.k',num2str(k),'.perc=0;']);
            eval(['S.Rs',num2str(source),'.k',num2str(k),'.pol=0;']);
            eval(['S.Rs',num2str(source),'.k',num2str(k),'.imp=1;']);
        end
    end
    eval(['S.Rs',num2str(source),'.amp=0;']);
    eval(['S.Rs',num2str(source),'.va=1;']);
    eval(['S.Rs',num2str(source),'.case.perc=100;']);
    eval(['S.Rs',num2str(source),'.case.pol=2;']);
end

% left sources
for source=1:4
    if contains(options.elmodel, 'Aleva')
        for k=12:23
            eval(['S.Ls',num2str(source),'.k',num2str(k),'.perc=0;']);
            eval(['S.Ls',num2str(source),'.k',num2str(k),'.pol=0;']);
            eval(['S.Ls',num2str(source),'.k',num2str(k),'.imp=1;']);
        end
    else
        for k=8:15
            eval(['S.Ls',num2str(source),'.k',num2str(k),'.perc=0;']);
            eval(['S.Ls',num2str(source),'.k',num2str(k),'.pol=0;']);
            eval(['S.Ls',num2str(source),'.k',num2str(k),'.imp=1;']);
        end
    end
    eval(['S.Ls',num2str(source),'.amp=0;']);
    eval(['S.Ls',num2str(source),'.va=1;']);
    eval(['S.Ls',num2str(source),'.case.perc=100;']);
    eval(['S.Ls',num2str(source),'.case.pol=2;']);
end

S.active=[1,1];
S.model='SimBio/FieldTrip (see Horn 2017)';
S.monopolarmodel=0;
S.amplitude={[0,0,0,0],[0,0,0,0]};
S=ea_activecontacts(S);
