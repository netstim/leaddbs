function [weight_r,rob_r,T,sig_T] = ea_stablecorr(y,X)

exponents=-2:0.5:2;
n=size(y,1);
allix=1:n;
sets=ea_kfoldsets(n,5,1);
% find most stable exponent for dataset
ea_dispercent(0,'Iterating exponents');
for exponent=1:length(exponents)
    for set=1:length(sets)
        othix=allix';
        othix(sets{set})=[];
        R(:,set)=ea_circcorr(y(othix),X(othix,:),exponents(exponent))';
    end
    cs=corr(R,'rows','pairwise');
    crosssims(exponent)=mean(cs(:));
    ea_dispercent(exponent/length(exponents));
end
ea_dispercent(1,'end');

[~,bestexp]=max(crosssims);
expn=exponents(bestexp);
disp(['Most stable exponent to apply = ',num2str(expn),'.']);

% now build stable mean:
clear R
for exponent=1:length(exponents)
    
    R(:,exponent)=ea_circcorr(y,X,exponents(exponent));
    
end

weights=crosssims.^5;
weights=weights./sum(weights);
weight_r=ea_nansum(R.*repmat(weights,size(R,1),1),2);

rob_r=ea_robustmean(R,2);

[~,p,~,stats]=ttest(atanh(R'));

T=stats.tstat';

sig_T=T;
sig_T(p>0.05)=nan;
