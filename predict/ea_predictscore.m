function [Ihat,Center]=ea_predictscore(XYZ,I,opts,cmd,cohs)
Center=nan; % center only returned if leave nothing out and distance as parameters.
if ~exist('opts','var')
    opts.usedist=1;
    opts.useexp='exp';
    opts.useregression='robust';
    opts.useweightedmean='lin';
    opts.discardnegativeweights=1;
end
warning off    
switch cmd
    case {'leaveOneOut','patient'}
        
        Ihat=ea_predictscore_leoo(XYZ,I,opts); % leave-one-out subroutine
        
    case {'leaveCohortOut','cohort'}
        Ihat=ea_predictscore_leco(XYZ,I,opts,cohs); % leave-one-cohort-out subroutine
        
    case {'leaveNothingOut','nothing'}
        [Ihat,Center]=dopredictscore(XYZ,I,1:size(XYZ,1),1:size(XYZ,1),opts);
end
warning on

function Ihat=ea_predictscore_leoo(XYZ,I,opts)
Ihat=nan(size(I));
N=size(XYZ,1);
for pt=1:N
    otherpts=1:N; otherpts(pt)=[];
    Ihat(pt)=dopredictscore(XYZ,I,otherpts,pt,opts);
end

function Ihat=ea_predictscore_leco(XYZ,I,opts,cohs)
Ihat=nan(size(I));
N=size(XYZ,1);
for coh=1:length(cohs)
    thiscoh=cohs{coh};
    othercohs=cohs; othercohs(coh)=[];
    otherpts=logical(sum(cell2mat(othercohs),2));
    Ihat(thiscoh)=dopredictscore(XYZ,I,otherpts,thiscoh,opts);
end


function [Ihat,Center]=dopredictscore(XYZ,I,modelpts,predictpts,opts)
Center=nan;

D=pdist2(XYZ,XYZ);
switch opts.useexp
    case 'exp'
        D=1./exp(D);
    case 'neg'
        D=-D;        
    case 'inv'
        D=1./D;
end


switch opts.useregression
    case 'svm'        
        Mdl = fitrsvm(XYZ((modelpts),:),I((modelpts),:),'Standardize',true,'KernelFunction','gaussian');
        Ihat = predict(Mdl,(XYZ(predictpts,:)));
    case {'distance','dist'}
        switch opts.useweightedmean 
            case 'off'
                Center=ea_nanmean(XYZ(modelpts,:));
            otherwise
                weights=getweights(I,modelpts,opts);
                Center=ea_nansum(XYZ(modelpts,:).*repmat(weights,1,3));
        end
        D=pdist2(XYZ(predictpts,:),Center);
        switch opts.useexp
            case 'exp'
                Ihat=1./exp(D);
            case 'neg'
                Ihat=-D;
            case 'inv'
                Ihat=1./D;
        end
    case 'regbeta'
        cnt=1; % we DO NEED to solve a GLM for each single patient!
        for pt=predictpts
            [b]=glmfit(D(modelpts,pt),I(modelpts));
            Ihat(cnt,1)=b(2);
            cnt=cnt+1;
        end
    case 'regR2'
        cnt=1; % we DO NEED to solve a GLM for each single patient!
        for pt=predictpts
            [b,bint,r,rint,stats]=regress(I(modelpts),ea_addone(D(modelpts,pt)));
            Ihat(cnt,1)=stats(1);
            cnt=cnt+1;
        end
    case 'regression'
    
        [b]=glmfit(D(modelpts,modelpts),I(modelpts));
        Ihat=ea_addone(D(predictpts,modelpts))*b;
        
    case 'regT'
        cnt=1; % we DO NEED to solve a GLM for each single patient!
        for pt=predictpts
            [b,dev,stats]=glmfit(D(modelpts,pt),I(modelpts));
            Ihat(cnt,1)=stats.t(2);
            cnt=cnt+1;
        end
    case 'robustbeta'
        cnt=1; % we DO NEED to solve a GLM for each single patient!
        for pt=predictpts
            [b]=robustfit(D(modelpts,pt),I(modelpts));
            Ihat(cnt,1)=b(2);
            cnt=cnt+1;
        end
    case 'robust'
        cnt=1; % we DO NEED to solve a GLM for each single patient!
        for pt=predictpts
            [b]=robustfit(D(modelpts,pt),I(modelpts));
            Ihat(cnt,1)=ea_addone(D(pt,pt))*b;
            cnt=cnt+1;
        end
        
    case 'robustcorr'
        cnt=1; % we DO NEED to solve a GLM for each single patient!
        for pt=predictpts
            [b,stats]=robustfit(D(modelpts,pt),I(modelpts));
            Ihat(cnt,1)=stats.coeffcorr(2);
            cnt=cnt+1;
        end
    case 'robustT'
        cnt=1; % we DO NEED to solve a GLM for each single patient!
        for pt=predictpts
            [b,stats]=robustfit(D(modelpts,pt),I(modelpts));
            Ihat(cnt,1)=stats.t(2);
            cnt=cnt+1;
        end
    case 'myskippedcorr'
        
        Ihat=myrobustcorr(D(modelpts,predictpts),I(modelpts));
        
    case 'regcorr'
        cnt=1; % we DO NEED to solve a GLM for each single patient!
        for pt=predictpts
            [b]=glmfit(D(modelpts,pt),I(modelpts));
            
            Ehat=ea_addone(D(modelpts,pt))*b;
            Ihat(cnt,1)=corr(Ehat,I(modelpts),'rows','pairwise','type','spearman');
            cnt=cnt+1;
        end
    case 'corr'
        Ihat=corr(D(modelpts,predictpts),I(modelpts),'rows','pairwise');
    case 'spearman'
        Ihat=corr(D(modelpts,predictpts),I(modelpts),'rows','pairwise','type','spearman');
    case 'kendall'
        Ihat=corr(D(modelpts,predictpts),I(modelpts),'rows','pairwise','type','kendall');
    case 'kendallcoeff'
        Ihat=KendallCoef([D(modelpts,predictpts),I(modelpts)]);
        
        Ihat=corr(D(modelpts,predictpts),I(modelpts),'rows','pairwise','type','kendall');
    case 'spearman_scatter'
        Ihat=ea_robustmean(corr(repmat(D(modelpts,predictpts),1,1000)+randn(length(modelpts),1000)*0.01,I(modelpts),'rows','pairwise','type','spearman'));

end
if opts.usedist && ~strcmp(opts.useregression,'distance')
    switch opts.useweightedmean
        case 'off'
            weights=ones(length(I(modelpts)),3);
        otherwise
            weights=getweights(I,modelpts,opts);
            weights=repmat(weights,1,3);
    end
    Dtomean=pdist2(ea_nanmean(XYZ(modelpts,:).*weights),XYZ(predictpts,:))';
    %Dtomean=pdist2(ea_nanmean(acs(modelpts,:)),acs(predictpts,:))';
    if opts.useexp==2 % use exp in both
        Dtomean=1./exp(Dtomean);
    else
        Dtomean=1./Dtomean;
        
        %        Dtomean=-Dtomean;
        %        Dtomean=Dtomean-min(Dtomean);
    end
    try
        Ihat=Ihat.*Dtomean;
    catch
        try
            Ihat=Ihat.*Dtomean';
        catch
            keyboard
        end
    end
end

function weights=getweights(I,modelpts,opts)
weights=I(modelpts);
if opts.discardnegativeweights
    weights(weights<0)=nan;
end
weights=weights-ea_nanmin(weights);
switch opts.useweightedmean
    case 'exp'
        weights=exp(weights);
    case 'square'
        weights=weights.^2;
end
weights=weights./ea_nansum(weights);


function fR=myrobustcorr(X,Y)

N=size(Y,1);
clear R
for leo=1:N
    all=1:N;
    all(rand(N,1)>0.5)=[];
    R(leo,:)=corr(X(all,:),Y(all),'rows','pairwise');
end
%outlsmag=R(logical(eye(N)'));
fR=ea_robustmean(R,1);
%figure, imagesc([corr(X,Y)';ea_robustmean(R,1)])

%Rs=corr([corr(X,Y)';ea_robustmean(R)]')

