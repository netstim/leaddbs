load('/path/to/your/LEAD_groupanalysis.mat'); % this creates variable M
M.ui.connectomename='HCP_MGH_30fold_groupconnectome (Horn 2017)'; % ... which can be modified afterwards
M.ui.listselect=1:length(M.patient.list); % usually helpful to create the model for all patients but then use only the ones defined in the script to run the analyses.
prefs=ea_prefs;

%% static part - usually no need to edit - can edit it by configuring lead group file correctly and closing it.
discfiberssetting = prefs.machine.lg.discfibers;
discfiberssetting.statmetric=1;

[fibsweighted,fibsin,fibsval,iaix]=ea_discfibers_calcdiscfibers(M,discfiberssetting);

%% flexible part - how to set up prediction - example is leave one cohort out crossvalidation:
Nperm=5000; % run as many as Nperm permutations
allpts=1:length(M.patient.list);
I=M.clinical.vars{M.ui.clinicallist};
nfibsval=fibsval; nfibsval(nfibsval==0)=nan; % only used in spearmans correlations

%nfibsval(nfibsval<50)=nan; % could be used to set a threshold on the E-Field method.
if discfiberssetting.statmetric==1
    % In t-test method, make sure fibers are at least connected to 0.2 percent of VTAs and not
    % more than to 0.8 percent.
    discthresh=0.2; % could be changed but should not exceed ~0.45
    fibsval(sum(fibsval,2)<round(discthresh*size(fibsval,2)),:)=0;
    fibsval(sum(fibsval,2)>round((1-discthresh)*size(fibsval,2)),:)=0;
end

Ihat = zeros(Nperm+1,length(I));
R0 = zeros(Nperm+1,1);

ea_dispercent(0,'Permutation testing');
for perm=1:Nperm+1
    switch discfiberssetting.statmetric
        case 1 % ttests, vtas - see Baldermann et al. 2019 Biological Psychiatry
            if perm>1
                Iperm=I(randperm(length(I)));
            else % use real empirical set in first run
                Iperm=I;
            end
            allvals=repmat(Iperm(allpts)',size(fibsval,1),1); % improvement values (taken from Lead group file or specified in line 12).
            
            fibsimpval=allvals; % Make a copy to denote improvements of connected fibers
            fibsimpval(~logical(fibsval))=nan; % Delete all unconnected values
            nfibsimpval=allvals; % Make a copy to denote improvements of unconnected fibers
            nfibsimpval(logical(fibsval))=nan; % Delete all connected values
            [~,p,~,Model]=ttest2(fibsimpval',nfibsimpval'); % Run two-sample t-test across connected / unconnected values
            Model.tstat(p>0.5)=nan; % discard noisy fibers (optional or could be adapted)
            for pt=allpts
                thisptval=fibsval(:,pt); % this patients connections to each fibertract (1 = connected, 0 = unconnected) 
                Ihat(perm,pt)=ea_nansum(Model.tstat'.*thisptval); % I hat is the estimate of improvements (not scaled to real improvements)
            end
        case 2 % spearmans correlations, efields - see Irmen et al. 2019 TBP
            Model=corr(nfibsval(:,allpts)',I(allpts),'rows','pairwise','type','Spearman'); % generate optimality values on all but left out patients
            for pt=allpts
                Ihat(perm,pt)=ea_nansum(Model.*nfibsval(:,pt)); % I hat is the estimate of improvements (not scaled to real improvements)
            end
    end
    
    R0(perm)=corr(Iperm,Ihat(perm,:)','type','Spearman','rows','pairwise');
    ea_dispercent(perm/Nperm);
end
ea_dispercent(1,'end');

R1=R0(1); % real correlation value when using empirical values
R0=R0(2:end); % 1-by-Nperm set of R values

% generate null distribution
R0=abs(R0);
R0=sort(R0,'descend');
p95=R0(round(0.05*Nperm));
v=ea_searchclosest(R0,R1);

pperm=v/Nperm;
disp(['Permuted p = ',sprintf('%0.2f',pperm),'.']);

h1=figure;
histogram(R0,round(Nperm/5),'FaceColor',[0.1,0.4,0.6],'EdgeColor',[0.1,0.4,0.6],'EdgeAlpha',0);
hold on
plot([R1,R1],[0,20],'Color',[0.6,0.2,0.3]);
text(R1,10,{'Unpermuted prediction \rightarrow',['p = ',sprintf('%0.2f',pperm)]},'Color',[0.6,0.2,0.3],'HorizontalAlignment','right');
saveas(h1,'my_result_permutation_test.png');



h2=ea_corrplot(I,Ihat(1,:)',{'Disc. Fiber prediction LOOCV','Empirical','Predicted'},'Spearman',M.patient.group);
saveas(h2,'my_result.png');






