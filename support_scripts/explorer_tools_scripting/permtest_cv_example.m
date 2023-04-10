clear
lead path
load('path/to/fibfilt.fibfilt','-mat'); % change to your fibfilt name

nperm=1000; % number of permutations
corrtype='Spearman'; % correlation type


% k-fold
for kfold=[5,10] % run it for 5 fold and 10-fold crossvalidations
    tractset.kfold=kfold;

    tractset.multitractmode = 'Single Tract Analysis';
    rng(tractset.rngseed);
    cvp = cvpartition(length(tractset.patientselection),'KFold',tractset.kfold);
    ea_dispercent(0,[num2str(kfold),'-fold: Iterating permutations']);
    R=zeros(nperm+1,1);
    for perm=0:nperm
        if perm % actual permutation run
            tractset.responsevar=tractset.responsevar(randperm(length(tractset.responsevar))); % permute improvements
            [I,Ihat] = crossval(tractset,cvp);
            R(perm+1)=corr(I,Ihat,'rows','pairwise','type',corrtype);
            ea_dispercent(0/nperm);
        else
            originalI=tractset.responsevar;
            [I,Ihat] = crossval(tractset,cvp);
            R(perm+1)=corr(I,Ihat,'rows','pairwise','type',corrtype);
        end
    end
    tractset.responsevar=originalI; % restore unpermuted improvements.
    ea_dispercent(1,'end');
    [h,p]=ea_plothistperm([num2str(kfold),'-fold: Null-Distribution'],R,{'Unpermuted Set'},{1},1,1);
    h.Position(3:4)=[1400,200]; % make figure long
    saveas(h,[num2str(kfold),'fold_nulldistribution.png']);
    saveas(h,[num2str(kfold),'fold_nulldistribution.fig']);
    saveas(h,[num2str(kfold),'fold_nulldistribution.eps']);
    close(h);
end

% leave one out

rng(tractset.rngseed);
cvp = cvpartition(length(tractset.patientselection), 'LeaveOut');
ea_dispercent(0,'Leave-One-Out: Iterating permutations');
R=zeros(nperm+1,1);
for perm=0:nperm
    if perm % actual permutation run
        tractset.responsevar=tractset.responsevar(randperm(length(tractset.responsevar))); % permute improvements
        [I,Ihat] = crossval(tractset,cvp);
        R(perm+1)=corr(I,Ihat,'rows','pairwise','type',corrtype);
        ea_dispercent(0/nperm);
    else
        originalI=tractset.responsevar;
        [I,Ihat] = crossval(tractset,cvp);
        R(perm+1)=corr(I,Ihat,'rows','pairwise','type',corrtype);
    end
end
tractset.responsevar=originalI; % restore unpermuted improvements.
ea_dispercent(1,'end');
[h,p]=ea_plothistperm([num2str(kfold),'-fold: Null-Distribution'],R,{'Unpermuted Set'},{1},1,1);
h.Position(3:4)=[1400,200]; % make figure long
saveas(h,'loo_nulldistribution.png');
saveas(h,'loo_nulldistribution.fig');
saveas(h,'loo_nulldistribution.eps');
close(h);