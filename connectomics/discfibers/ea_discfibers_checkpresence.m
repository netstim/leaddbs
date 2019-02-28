function [reforce,connectomechanged,reformat]=ea_discfibers_checkpresence(M,opts)
reforce=1; connectomechanged=1; reformat=1;

switch opts.statmetric
    case 1 % ttests
        savesuffix='_ttests';
    case 2 % spearmans R
        savesuffix='_spearmansrho';
end

if M.ui.mirrorsides
    msuffix='_mirrored';
else
    msuffix='';
end

if exist([M.ui.groupdir,'correlative_fibertracts',msuffix,savesuffix,'.mat'],'file')
    d=load([M.ui.groupdir,'correlative_fibertracts',msuffix,savesuffix,'.mat'],'opts');
    if isequaln(opts,d.opts)
        reforce=0;
    end
end

if exist([M.ui.groupdir,'connected_fibers',msuffix,'.mat'],'file') % check if base connectome changed.
    d=load([M.ui.groupdir,'connected_fibers',msuffix,'.mat'],'opts');
    if isequaln(d.opts.connectome,opts.connectome)
        connectomechanged=0;
    end
end

if ~reforce
    if exist([M.ui.groupdir,'correlative_fibertracts_reformatted',msuffix,savesuffix,'.mat'],'file') % check if base connectome changed.
        d=load([M.ui.groupdir,'correlative_fibertracts_reformatted',msuffix,savesuffix,'.mat'],'opts');
        if isequaln(d.opts.connectome,opts.connectome) && isequaln(d.opts.statmetric,opts.statmetric)
            reformat=0;
        end
    end
end
