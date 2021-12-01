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

for p=1:length(PL)
    if isfield(PL(p),'vatsurfs')
        delete(PL(p).vatsurfs);
    end

    if isfield(PL(p),'matsurf')
        delete(PL(p).matsurf);
    end

    if isfield(PL(p),'matseedsurf')
        delete(PL(p).matseedsurf);
    end

    if isfield(PL(p),'matsurf')
        delete(PL(p).matsurf);
    end

    if isfield(PL(p),'quiv')
        delete(PL(p).quiv);
    end

    if isfield(PL(p),'fib_plots')
        if isfield(PL(p).fib_plots,'fibs')
            delete(PL(p).fib_plots.fibs);
        end

        if isfield(PL(p).fib_plots,'dcfibs')
            delete(PL(p).fib_plots.dcfibs);
        end
    end

    if isfield(PL(p),'regionsurfs')
        delete(PL(p).regionsurfs);
    end

    if isfield(PL(p),'conlabels')
        delete(PL(p).conlabels);
    end

    if isfield(PL(p),'ht')
        delete(PL(p).ht);
    end

    if isfield(PL(p),'fiberActivation')
        cellfun(@delete, PL(p).fiberActivation);
    end
end
