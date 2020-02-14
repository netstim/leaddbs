function ea_deletePL(resultfig,varname,mode)


PL=getappdata(resultfig,[mode,varname]);
setappdata(resultfig,[mode,'fibersfile'],'nan');
setappdata(resultfig,[mode,'seedfile'],'nan');
setappdata(resultfig,[mode,'targetsfile'],'nan');
setappdata(resultfig,[mode,'thresh'],'nan');

try
if verLessThan('matlab','8.5') % ML <2014a support
    for p=1:length(PL)
        if isfield(PL(p),'vatsurfs')
            delete(PL(p).vatsurfs(logical(PL(p).vatsurfs)));
        end
        if isfield(PL(p),'quiv')
            delete(PL(p).quiv(logical(PL(p).quiv)));
        end   
        if isfield(PL(p),'matseedsurf')
            for t=1:length(PL(p).matseedsurf)
                try delete(PL(p).matseedsurf{t});
                end
            end
        end
        if isfield(PL(p),'matsurf')
            for t=1:length(PL(p).matsurf)
                try delete(PL(p).matsurf{t}); end
            end
        end
        if isfield(PL(p),'fib_plots')
            if isfield(PL(p).fib_plots,'fibs')
                delete(PL(p).fib_plots.fibs(logical(PL(p).fib_plots.fibs)));
            end
            
            if isfield(PL(p).fib_plots,'dcfibs')
                todelete=PL(p).fib_plots.dcfibs((PL(p).fib_plots.dcfibs(:)>0));
                delete(todelete(:));
                
            end
        end
        if isfield(PL(p),'regionsurfs')
            todelete=PL(p).regionsurfs(logical(PL(p).regionsurfs));
            delete(todelete(:));
        end
        
        if isfield(PL(p),'conlabels')
            todelete=PL(p).conlabels(logical(PL(p).conlabels));
            delete(todelete(:));
        end
        if isfield(PL(p),'ht')
            delete(PL(p).ht);
        end
    end
    
    
else
    for p=1:length(PL) 
        
        if isfield(PL(p),'matseedsurf')
            for t=1:length(PL(p).matseedsurf)
                try delete(PL(p).matseedsurf{t});
                end
            end
        end
        if isfield(PL(p),'matsurf')
            for t=1:length(PL(p).matsurf)
                try delete(PL(p).matsurf{t}); end
            end
        end
        if isfield(PL(p),'vatsurfs')
            delete(PL(p).vatsurfs);
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
    end
    
end
end