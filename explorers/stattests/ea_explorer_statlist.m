function tests = ea_explorer_statlist(varargin)

% first input can be the variable, if applied, will only return the correct
% tests (e.g. binary / continuous).

if nargin
    outcomein=varargin{1};
    if iscell(outcomein)
        outcomein=outcomein{1};
    end
    isbinary=ea_isbinary(outcomein);
else
    isbinary=nan;
end

    warning('off','all')
    tests = table('Size',[0,5],...
        'VariableTypes',{'categorical','string','categorical','cell','cell'},...
        'VariableNames',{'name','file','type','outcometype','compatibility'});
    row=0;

[pth]=fileparts(mfilename('fullpath')); % get dir
funs=dir(fullfile(pth,'ea_explorer_stats_*.m'));
cnt=1;
for fun=1:length(funs)
entry=feval(ea_stripext(funs(fun).name));

if isnan(isbinary)
    tests.name(cnt)=entry.name;
    tests.file(cnt)=entry.file;
    tests.type(cnt)=entry.type;
    tests.outcometype{cnt}=entry.outcometype;
    tests.compatibility{cnt}=entry.compatibility;
    cnt=cnt+1;
else
    if isbinary
        if ismember('binary',entry.outcometype) % only add test to list if supported for binary variables
            tests.name(cnt)=entry.name;
            tests.file(cnt)=entry.file;
            tests.type(cnt)=entry.type;
            tests.outcometype{cnt}=entry.outcometype;
            tests.compatibility{cnt}=entry.compatibility;
            cnt=cnt+1;
        end
    end
    if ~isbinary
        if ismember('gradual',entry.outcometype) % only add test to list if supported for binary variables
            tests.name(cnt)=entry.name;
            tests.file(cnt)=entry.file;
            tests.type(cnt)=entry.type;
            tests.outcometype{cnt}=entry.outcometype;
            tests.compatibility{cnt}=entry.compatibility;
            cnt=cnt+1;
        end
    end
end

end