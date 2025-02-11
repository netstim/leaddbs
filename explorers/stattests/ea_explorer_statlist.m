function tests = ea_explorer_statlist(varargin)
% first input can be the variable, if applied, will only return the correct
% tests (e.g. binary / continuous).

if nargin
    outcomein = varargin{1};
    if iscell(outcomein)
        outcomein = outcomein{1};
    end
    outcomein(isnan(outcomein)) = 0; %convert NaN in your outcomes to 0 because logical cannot handle NaN values
    isbinary = ea_isbinary(outcomein);
else
    isbinary = nan;
end

oldWarningState = warning;
warning('off', 'all')

tests = table('Size',[0,5],...
    'VariableTypes',{'categorical','string','categorical','cell','cell'},...
    'VariableNames',{'name','file','type','outcometype','compatibility'});

funcs = dir(fullfile(fileparts(mfilename('fullpath')), 'ea_explorer_stats_*.m'));
funcs = ea_stripext({funcs.name}');

cnt = 1;
for i=1:length(funcs)
    entry = feval(funcs{i});
    if isnan(isbinary)
        tests.name(cnt) = entry.name;
        tests.file(cnt) = entry.file;
        tests.type(cnt) = entry.type;
        tests.outcometype{cnt} = entry.outcometype;
        tests.compatibility{cnt} = entry.compatibility;
        cnt = cnt+1;
    else
        if isbinary
            if ismember('binary', entry.outcometype) % only add test to list if supported for binary variables
                tests.name(cnt) = entry.name;
                tests.file(cnt) = entry.file;
                tests.type(cnt) = entry.type;
                tests.outcometype{cnt} = entry.outcometype;
                tests.compatibility{cnt} = entry.compatibility;
                cnt = cnt+1;
            end
        end
        if ~isbinary
            if ismember('gradual', entry.outcometype) % only add test to list if supported for binary variables
                tests.name(cnt) = entry.name;
                tests.file(cnt) = entry.file;
                tests.type(cnt) = entry.type;
                tests.outcometype{cnt} = entry.outcometype;
                tests.compatibility{cnt} = entry.compatibility;
                cnt = cnt+1;
            end
        end
    end
end

warning(oldWarningState);
