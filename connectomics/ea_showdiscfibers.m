function ea_showdiscfibers(M,discfiberssetting,resultfig)

patlist=M.patient.list(M.ui.listselect);
I=M.clinical.vars{M.ui.clinicallist}(M.ui.listselect);

thresh=150;% should only apply for heatfibertracts_corr
tic

% Get discriminative fiber setting
connthreshold = discfiberssetting.connthreshold/100;
showfibersset = discfiberssetting.showfibersset;
pospredthreshold = discfiberssetting.pospredthreshold/100;
negpredthreshold = discfiberssetting.negpredthreshold/100;
statmetric=discfiberssetting.statmetric;

% protocol selection to be able to check if same analysis has been run
% before.
%opts.percent=connthreshold; % need not protocol anymore, is now dynamic
%option.
opts.patientselection=M.ui.listselect;
opts.regressor=I;
opts.connectome=M.ui.connectomename;
opts.allpatients=M.patient.list;
opts.mirrorsides=M.ui.mirrorsides;
opts.statmetric=statmetric;

if M.ui.mirrorsides
    msuffix='_mirrored';
else
    msuffix='';
end

[reforce,connectomechanged,reformat]=checkpresence(M,opts); % only static opts need to be equal.
if reforce
    allroilist=cell(length(M.patient.list),2);
    switch statmetric
        case 1 % use paired T-Tests and binary VTA
            suffix='';
        case 2 % use Spearman Rs and E-Fields
            suffix='_efield';
            prefs=ea_prefs;
            if strcmp(prefs.lcm.vatseed,'efield_gauss')
                suffix='_efield_gauss';
            end
    end
    
    for sub=1:length(M.patient.list) % all patients - for connected fibers selection
        if ~M.ui.mirrorsides
            allroilist{sub,1}=[M.patient.list{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat',suffix,'_right.nii'];
            allroilist{sub,2}=[M.patient.list{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat',suffix,'_left.nii'];
        else
            allroilist(sub,:)=ea_genflippedjointnii([M.patient.list{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat',suffix,'_right.nii'],[M.patient.list{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat',suffix,'_left.nii'],statmetric==1);
        end
    end

    if ~exist([M.ui.groupdir,'connected_fibers',msuffix,'.mat'],'file')
        cfile=[ea_getconnectomebase('dMRI'),M.ui.connectomename,filesep,'data.mat'];
    else
        if connectomechanged
            cfile=[ea_getconnectomebase('dMRI'),M.ui.connectomename,filesep,'data.mat'];
        else
            cfile=[M.ui.groupdir,'connected_fibers',msuffix,'.mat'];
        end
    end
    switch statmetric
        case 1 % ttests
            [fibsweighted,fibsin]=ea_heatfibertracts(cfile,{allroilist},M.ui.listselect,{I},thresh,connthreshold);
        case 2 % spearmans R
            [fibsweighted,fibsin]=ea_heatfibertracts_corr(cfile,{allroilist},M.ui.listselect,{I},thresh,connthreshold);
    end
    save([M.ui.groupdir,'connected_fibers',msuffix,'.mat'],'fibsin','opts','-v7.3');
    save([M.ui.groupdir,'correlative_fibertracts',msuffix,'.mat'],'fibsweighted','opts','-v7.3');
else
    load([M.ui.groupdir,'correlative_fibertracts',msuffix,'.mat']);
end

% visualize:
if reformat
    fibidx=unique(fibsweighted(:,4));
    fibcell=cell(length(fibidx),1);
    vals=zeros(length(fibidx),1);
    cnt=1;
    ea_dispercent(0,'reformatting fibers');
    for fib=fibidx'
        fibcell{cnt}=fibsweighted(fibsweighted(:,4)==fib,1:3);
        vals(cnt)=mean(fibsweighted(fibsweighted(:,4)==fib,5));
        cnt=cnt+1;
        ea_dispercent(cnt/length(fibidx));
    end
    ea_dispercent(1,'end');
    save([M.ui.groupdir,'correlative_fibertracts_reformatted',msuffix,'.mat'],'fibcell','vals','opts','-v7.3');
else
    load([M.ui.groupdir,'correlative_fibertracts_reformatted',msuffix,'.mat']);
end

set(0,'CurrentFigure',resultfig);

% Normalize vals
vals(isnan(vals))=0;
vals=vals./max(abs(vals));

% vals and fibcell to be trimmed for visualization
tvals=vals;
tfibcell=fibcell;

% Calculate positive/negative threshold for positive/negative predictive
% fibers according to 'predthreshold'
posits=tvals(tvals>0);
negits=tvals(tvals<0);
posits=sort(posits,'descend');
negits=sort(negits,'ascend');
posthresh=posits(round(length(posits)*pospredthreshold));
negthresh=negits(round(length(negits)*negpredthreshold));

% Save the original values for reusing in slider
setappdata(resultfig, 'vals', vals);
setappdata(resultfig, 'fibcell', fibcell);
setappdata(resultfig, 'showfibersset', showfibersset);
setappdata(resultfig, 'pospredthreshold', pospredthreshold);
setappdata(resultfig, 'negpredthreshold', negpredthreshold);
setappdata(resultfig, 'posits', posits);
setappdata(resultfig, 'negits', negits);

switch showfibersset
    case 'positive'
        negthresh = negits(1)-eps;
        disp(['Fiber colors: Positive (T = ',num2str(posthresh),' ~ ',num2str(posits(1)), ')']);
    case 'negative'
        posthresh = posits(1)+eps;
        disp(['Fiber colors: Negative (T = ',num2str(negits(1)),' ~ ',num2str(negthresh), ')']);
    case 'both'
        disp(['Fiber colors: Positive (T = ',num2str(posthresh),' ~ ',num2str(posits(1)), ...
          '); Negative (T = ',num2str(negits(1)),' ~ ',num2str(negthresh),').']);
end

% Remove tvals and fibers outside the thresholding range
remove=logical(logical(tvals<posthresh) .* logical(tvals>negthresh));
tvals(remove)=[];
tfibcell(remove)=[];

% Rescale positive/negative tvals to [0 1]/[-1 0]
tvalsRescale = tvals;
tvalsRescale(tvals>0)=ea_minmax(tvals(tvals>0));
tvalsRescale(tvals<0)=-ea_minmax(-tvals(tvals<0));

% Contruct colormap
colormap gray
map=ea_redblue(1024);
fibcolorInd=tvalsRescale*(size(map,1)/2-0.5);
fibcolorInd=fibcolorInd+(size(map,1)/2+0.5);

% Set alphas of fibers with light color to 0
colorbarThreshold = 0.60; % Percentage of the pos/neg color to be kept
negUpperBound=ceil(size(map,1)/2*colorbarThreshold);
poslowerBound=floor((size(map,1)-size(map,1)/2*colorbarThreshold));
alphas=zeros(size(fibcolorInd,1),1);
switch showfibersset
    case 'positive'
        alphas(round(fibcolorInd)>=poslowerBound) = 1;
    case 'negative'
        alphas(round(fibcolorInd)<=negUpperBound) = 1;
    case 'both'
        alphas(round(fibcolorInd)>=poslowerBound) = 1;
        alphas(round(fibcolorInd)<=negUpperBound) = 1;
end

alphas(round(fibcolorInd)>=poslowerBound) = 1;
fibalpha=mat2cell(alphas,ones(size(fibcolorInd,1),1));

% Plot fibers
h=streamtube(tfibcell,0.2);
nones=repmat({'none'},size(fibcolorInd));
[h.EdgeColor]=nones{:};

% Calulate fiber colors
colors=map(round(fibcolorInd),:);
fibcolor=mat2cell(colors,ones(size(fibcolorInd)));

% Set fiber colors and alphas
[h.FaceColor]=fibcolor{:};
[h.FaceAlpha]=fibalpha{:};

setappdata(resultfig, 'discfibers', h);

% Set colorbar tick positions and labels
cbvals = tvals(logical(alphas));
% cbvals=tvalsRescale(logical(alphas));
switch showfibersset
    case 'positive'
        cbmap = map(ceil(length(map)/2+0.5):end,:);
        tick = [poslowerBound, length(map)] - floor(length(map)/2) ;
        poscbvals = sort(cbvals(cbvals>0));
        ticklabel = [poscbvals(1), poscbvals(end)];
        ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
    case 'negative'
        cbmap = map(1:floor(length(map)/2-0.5),:);
        tick = [1, negUpperBound];
        negcbvals = sort(cbvals(cbvals<0));
        ticklabel = [negcbvals(1), negcbvals(end)];
        ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
    case 'both'
        cbmap = map;
        tick = [1, negUpperBound, poslowerBound, length(map)];
        poscbvals = sort(cbvals(cbvals>0));
        negcbvals = sort(cbvals(cbvals<0));
        ticklabel = [min(cbvals), negcbvals(end), poscbvals(1), max(cbvals)];
        ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
end

% Plot colorbar
cbfig = ea_plot_colorbar(cbmap, [], 'v', '', tick, ticklabel);
set(cbfig, 'NumberTitle', 'off');
setappdata(resultfig, 'cbfig', cbfig);

% Discriminative fiber control
discfiberscontrol = ea_discfiberscontrol(resultfig);
setappdata(resultfig, 'discfiberscontrol', discfiberscontrol);





function [reforce,connectomechanged,reformat]=checkpresence(M,opts)
reforce=1; connectomechanged=1; reformat=1;
if M.ui.mirrorsides
    msuffix='_mirrored';
else
    msuffix='';
end

if exist([M.ui.groupdir,'correlative_fibertracts',msuffix,'.mat'],'file')
    d=load([M.ui.groupdir,'correlative_fibertracts',msuffix,'.mat'],'opts');
    if isequaln(opts,d.opts)
        reforce=0;
    end
end

if exist([M.ui.groupdir,'connected_fibers',msuffix,'.mat'],'file') % check if base connectome changed.
    d=load([M.ui.groupdir,'correlative_fibertracts',msuffix,'.mat'],'opts');
    if isequaln(d.opts.connectome,opts.connectome)
        connectomechanged=0;
    end
end

if ~reforce
    if exist([M.ui.groupdir,'correlative_fibertracts_reformatted',msuffix,'.mat'],'file') % check if base connectome changed.
        d=load([M.ui.groupdir,'correlative_fibertracts_reformatted',msuffix,'.mat'],'opts');
        if isequaln(d.opts.connectome,opts.connectome) && isequaln(d.opts.statmetric,opts.statmetric)
            reformat=0;
        end
    end
end


function str=pointtodash(str)
str=strrep(str,'.','-');


function str=stripblanks(str)
str=strrep(str,'(','');
str=strrep(str,')','');
str=strrep(str,' ','');
