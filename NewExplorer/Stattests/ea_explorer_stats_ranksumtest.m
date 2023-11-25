function [valsout,psout]=ea_explorer_stats_ranksumtest(valsin,outcomein)
outcomein=repmat(outcomein',size(valsin,1),1);

group1=valsin .* outcomein;
group2=double(isnan(valsin));
group2(group2==0)=nan;
group2=group2 .* outcomein;

psout=nan(size(valsin,1),1);
for i=1:size(valsin,1)
    [psout(i),~,stats(i)]=ranksum(group1(i,:),group2(i,:));
end
zval=[stats(:).zval]';
valsout = sign(zval).*-log10(psout);
end