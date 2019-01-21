function ea_showdiscfibers(M,discfiberssetting,resultfig)

patlist=M.patient.list(M.ui.listselect);
I=M.clinical.vars{M.ui.clinicallist}(M.ui.listselect);

thresh=0.5;
tic

% Get discriminative fiber setting
connthreshold = discfiberssetting.connthreshold/100;
predthreshold = discfiberssetting.predthreshold/100;
showpositiveonly = discfiberssetting.showpositiveonly;

% protocol selection to be able to check if same analysis has been run
% before.
opts.percent=connthreshold;
opts.patientselection=M.ui.listselect;
opts.regressor=I;
opts.connectome=M.ui.connectomename;
opts.allpatients=M.patient.list;

[reforce,connectomechanged,reformat]=checkpresence(M,opts);
if reforce
    allroilist=cell(length(M.patient.list),2);
    for sub=1:length(M.patient.list) % all patients - for connected fibers selection
        allroilist{sub,1}=[M.patient.list{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat_right.nii'];
        allroilist{sub,2}=[M.patient.list{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat_left.nii'];
    end

    if ~exist([M.ui.groupdir,'connected_fibers.mat'],'file')
        cfile=[ea_getconnectomebase('dMRI'),M.ui.connectomename,filesep,'data.mat'];
    else
        if connectomechanged
            cfile=[ea_getconnectomebase('dMRI'),M.ui.connectomename,filesep,'data.mat'];
        else
            cfile=[M.ui.groupdir,'connected_fibers.mat'];
        end
    end
    [fibsweighted,fibsin]=ea_heatfibertracts(cfile,{allroilist},M.ui.listselect,{I},thresh,connthreshold);
    save([M.ui.groupdir,'connected_fibers.mat'],'fibsin','opts','-v7.3');
    save([M.ui.groupdir,'correlative_fibertracts.mat'],'fibsweighted','opts','-v7.3');
else
    load([M.ui.groupdir,'correlative_fibertracts.mat']);
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
    save([M.ui.groupdir,'correlative_fibertracts_reformatted.mat'],'fibcell','vals','opts','-v7.3');
else
    load([M.ui.groupdir,'correlative_fibertracts_reformatted.mat']);
end

set(0,'CurrentFigure',resultfig);

% tvals for thresholding
tvals=vals;
tvals(isnan(tvals))=0;
tvals=tvals./max(abs(tvals));

% Calculate positive/negative threshold for positive/negative predictive
% fibers according to 'predthreshold'
posits=tvals(tvals>0);
negits=tvals(tvals<0);
posits=sort(posits,'descend');
negits=sort(negits,'ascend');
posthresh=posits(round(length(posits)*predthreshold));

if showpositiveonly
    negthresh = negits(1)-eps;
    disp(['Fiber colors: Positive (T = ',num2str(posthresh),' ~ ',num2str(max(posits)), ')']);
else
    negthresh=negits(round(length(negits)*predthreshold));
    disp(['Fiber colors: Positive (T = ',num2str(posthresh),' ~ ',num2str(max(posits)), ...
          '); Negative (T = ',num2str(min(negits)),' ~ ',num2str(negthresh),').']);
end

% Remove tvals and fibers outside the thresholding range
remove=logical(logical(tvals<posthresh) .* logical(tvals>negthresh));
tvals(remove)=[];
fibcell(remove)=[];

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
if ~showpositiveonly
    alphas(round(fibcolorInd)<=negUpperBound) = 1;
end
alphas(round(fibcolorInd)>=poslowerBound) = 1;
fibalpha=mat2cell(alphas,ones(size(fibcolorInd,1),1));

% Plot fibers
h=streamtube(fibcell,0.2);
nones=repmat({'none'},size(fibcolorInd));
[h.EdgeColor]=nones{:};

% Calulate fiber colors
colors=map(round(fibcolorInd),:);
fibcolor=mat2cell(colors,ones(size(fibcolorInd)));

% Set fiber colors and alphas
[h.FaceColor]=fibcolor{:};
[h.FaceAlpha]=fibalpha{:};

% Set colorbar tick positions and labels
cbvals = tvals(logical(alphas));
% cbvals=tvalsRescale(logical(alphas));
if showpositiveonly
    cbmap = map(ceil(length(map)/2+0.5):end,:);
    tick = [poslowerBound, length(map)] - floor(length(map)/2) ;
    poscbvals = sort(cbvals(cbvals>0));
    ticklabel = [poscbvals(1), max(cbvals)];
    ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
else
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


function [reforce,connectomechanged,reformat]=checkpresence(M,opts)
reforce=1; connectomechanged=1; reformat=1;
if exist([M.ui.groupdir,'correlative_fibertracts.mat'],'file')
    d=load([M.ui.groupdir,'correlative_fibertracts.mat'],'opts');
    if isequal(opts,d.opts)
        reforce=0;
    end
end

if exist([M.ui.groupdir,'connected_fibers.mat'],'file') % check if base connectome changed.
    d=load([M.ui.groupdir,'correlative_fibertracts.mat'],'opts');
    if isequal(d.opts.connectome,opts.connectome)
        connectomechanged=0;
    end
end

if ~reforce
    if exist([M.ui.groupdir,'correlative_fibertracts_reformatted.mat'],'file') % check if base connectome changed.
        d=load([M.ui.groupdir,'correlative_fibertracts_reformatted.mat'],'opts');
        if isequal(d.opts.connectome,opts.connectome)
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
