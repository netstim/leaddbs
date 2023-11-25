function [valsout,psout]=ea_explorer_stats_2samplettest(valsin,outcomein)
outcomein=repmat(outcomein',size(valsin,1),1);

group1=valsin .* outcomein;
group2=double(isnan(valsin));
group2(group2==0)=nan;
group2=group2 .* outcomein;
valsout=nan(size(valsin,1),1);
psout=nan(size(valsin,1),1);
[~,psout,~,stats]=ttest2(group1',group2');
psout=psout';
valsout = stats.tstat';
end