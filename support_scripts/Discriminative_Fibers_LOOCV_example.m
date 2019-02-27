load('/path/to/your/LEAD_groupanalysis.mat'); % this creates variable M
M.ui.connectomename='HCP_MGH_30fold_groupconnectome (Horn 2017)'; % ... which can be modified afterwards
prefs=ea_prefs;

%% static part - usually no need to edit - can edit it by configuring lead group file correctly and closing it.
discfiberssetting = prefs.machine.lg.discfibers;
[fibsweighted,fibsin,fibsval,iaix]=ea_calcdiscfibers(M,discfiberssetting);


%% flexible part - how to set up prediction - example is leave one out crossvalidation:
allpts=1:length(M.patient.list);
I=M.clinical.vars{M.ui.clinicallist};
ea_dispercent(0,'Predicting left-out patients');
for pt=allpts
    opts=allpts; opts(opts==pt)=[];
    switch discfiberssetting.statmetric
        case 1 % ttests, vtas
            thisptval=fibsval(:,pt);
            optsval=fibsval(:,opts);
            allvals=repmat(I',size(optsval,1),1);
            fibsimpval=allvals;
            fibsimpval(~logical(optsval))=nan;
            nfibsimpval=allvals;
            nfibsimpval(logical(optsval))=nan;
            [~,~,~,Model]=ttest2(fibsimpval',nfibsimpval');
            Ihat(pt)=ea_nansum(Model.tstat'.*thisptval); % I hat is the estimate of improvements (not scaled to real improvements)
        case 2 % spearmans correlations, efields
            Model=corr(fibsval(:,opts)',I(opts),'rows','pairwise','type','Spearman');
            Ihat(pt)=ea_nansum(Model.*nfibsval(:,pt)); % I hat is the estimate of improvements (not scaled to real improvements)
    end
    ea_dispercent(pt/length(allpts));
end
ea_dispercent(1,'end');

ea_corrplot(I,Ihat',{'Disc. Fiber prediction',I,Ihat},'permutation_spearman');







