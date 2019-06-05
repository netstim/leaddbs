function lc=ea_initlcopts(varargin)
% usage: ea_initlcopts([handletolcfig],[lc])

if nargin>1
    lc=varargin{2};
else
    if isempty(varargin{1})
        h=lead_connectome;
    else
        h=varargin{1};
    end
    lc.general.parcellation=getappdata(h,'parcellation');
    lc.general.parcellationn=1;
    lc.general.parcellation=lc.general.parcellation{lc.general.parcellationn};
    lc.graph.degree_centrality=0;
    lc.graph.eigenvector_centrality=0;
    lc.graph.nodal_efficiency=0;
    lc.graph.struc_func_sim=0;
    lc.graph.fthresh=nan;
    lc.graph.sthresh=nan;
    lc.func.compute_CM=0;
    lc.func.compute_GM=0;
    lc.func.prefs.TR=2;
    lc.struc.compute_CM=0;
    lc.struc.compute_GM=0;
    lc.struc.ft.do=0;
    ftmethods=getappdata(h,'ftmethod');
    if isempty(varargin{1})
        close(h)
    end
    lc.struc.ft.method='ea_ft_gqi_yeh';
    lc.struc.ft.methodn=find(ismember(ftmethods, 'ea_ft_gqi_yeh'));
    lc.struc.ft.dsistudio.fiber_count=200000;
    lc.struc.ft.normalize=0;
end
lc.nbs.adv.method=1;
lc.nbs.adv.compsize=1;
lc.nbs.adv.perm=5000;
lc.nbs.adv.alpha=0.05;
lc.nbs.adv.exch='';
