function ea_discfibers_showdiscfibers(M,discfiberssetting,resultfig,fibsweighted)

I=M.clinical.vars{M.ui.clinicallist}(M.ui.listselect);

% Get discriminative fiber setting
statmetric = discfiberssetting.statmetric;

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

patlist=M.patient.list(M.ui.listselect);
showfibersset = discfiberssetting.showfibersset;
pospredthreshold = discfiberssetting.pospredthreshold/100;
negpredthreshold = discfiberssetting.negpredthreshold/100;

[reforce,connectomechanged,reformat]=ea_discfibers_checkpresence(M,opts); % only static opts need to be equal.
switch statmetric
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
    save([M.ui.groupdir,'correlative_fibertracts_reformatted',msuffix,savesuffix,'.mat'],'fibcell','vals','opts','-v7.3');
else
    load([M.ui.groupdir,'correlative_fibertracts_reformatted',msuffix,savesuffix,'.mat']);
end

set(0,'CurrentFigure',resultfig);

% Normalize vals
vals(isnan(vals))=0;
% vals=vals./max(abs(vals));

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

discfiberID = '';

% Save the original values for reusing in slider
setappdata(resultfig, ['vals',discfiberID], vals);
setappdata(resultfig, ['fibcell',discfiberID], fibcell);
setappdata(resultfig, ['showfibersset',discfiberID], showfibersset);
setappdata(resultfig, ['pospredthreshold',discfiberID], pospredthreshold);
setappdata(resultfig, ['negpredthreshold',discfiberID], negpredthreshold);
setappdata(resultfig, ['posits',discfiberID], posits);
setappdata(resultfig, ['negits',discfiberID], negits);

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
tvalsRescale(tvals>0)=ea_rescale(tvals(tvals>0), [0 1]);
tvalsRescale(tvals<0)=ea_rescale(tvals(tvals<0), [-1 0]);

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

% only show right side
% for i=1:length(tfibcell)
%     if any(tfibcell{i}(:,1)<0)
%         alphas(i) = 0;
%     end
% end

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

figTitle = 'discfibers';
discfibersname = ['discfibers', discfiberID];
cbfigname = ['cbfig', discfiberID];
discfiberscontrolname = ['discfiberscontrol', discfiberID];

% Plot colorbar
cbfig = ea_plot_colorbar(cbmap, [], 'h', '', tick, ticklabel);
set(cbfig, 'NumberTitle', 'off',  'Name', ['Colorbar: ', figTitle]);

% Discriminative fiber control
discfiberscontrol = ea_discfibers_control(resultfig, discfiberID);
set(discfiberscontrol, 'NumberTitle', 'off', 'Name', ['Control: ', figTitle]);
setappdata(discfiberscontrol, 'discfiberID', discfiberID);

setappdata(resultfig, discfibersname, h);
setappdata(resultfig, cbfigname, cbfig);
setappdata(resultfig, discfiberscontrolname, discfiberscontrol);
set(0,'CurrentFigure',resultfig)


function str=pointtodash(str)
str=strrep(str,'.','-');


function str=stripblanks(str)
str=strrep(str,'(','');
str=strrep(str,')','');
str=strrep(str,' ','');
