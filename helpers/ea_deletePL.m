function ea_deletePL(varargin)

if nargin == 1 % Used in ea_stimparams.m
    PL = varargin{1};
else % Used in convis
    resultfig = varargin{1};
    varname = varargin{2};
    mode = varargin{3};
    PL = getappdata(resultfig, [mode,varname]);
    setappdata(resultfig, [mode,'fibersfile'], 'nan');
    setappdata(resultfig, [mode,'seedfile'], 'nan');
    setappdata(resultfig, [mode,'targetsfile'], 'nan');
    setappdata(resultfig, [mode,'thresh'], 'nan');
end
if ~isempty(PL)
    for p=1:length(PL)
        if isfield(PL(p),'vatsurfs')
            try delete(PL(p).vatsurfs); end
        end

        if isfield(PL(p),'matsurf')
            if iscell(PL(p).matsurf)
                for c=1:length(PL(p).matsurf)
                    try delete(PL(p).matsurf{c}); end
                end
            else
                try delete(PL(p).matsurf); end
            end
        end

        if isfield(PL(p),'matseedsurf')
            if iscell(PL(p).matseedsurf)
                for c=1:length(PL(p).matseedsurf)
                    try delete(PL(p).matseedsurf{c}); end
                end
            else
                try delete(PL(p).matseedsurf); end
            end
        end



        if isfield(PL(p),'quiv')
            try delete(PL(p).quiv); end
        end

        if isfield(PL(p),'fib_plots')
            if isfield(PL(p).fib_plots,'fibs')
                try delete(PL(p).fib_plots.fibs); end
            end

            if isfield(PL(p).fib_plots,'dcfibs')
                try delete(PL(p).fib_plots.dcfibs); end
            end
        end

        if isfield(PL(p),'regionsurfs')
            try delete(PL(p).regionsurfs); end
        end

        if isfield(PL(p),'conlabels')
            if iscell(PL(p).conlabels)
                for c=1:length(PL(p).conlabels)
                    try                    delete(PL(p).conlabels{c}); end
                end
            else
                try delete(PL(p).conlabels); end
            end
        end

        if isfield(PL(p),'ht')
            try delete(PL(p).ht); end
        end

        if isfield(PL(p),'fiberActivation')
            cellfun(@delete, PL(p).fiberActivation);
        end
    end
end
