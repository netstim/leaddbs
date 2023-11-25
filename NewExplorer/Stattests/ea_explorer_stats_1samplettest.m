function [valsout,psout]=ea_explorer_stats_1samplettest(valsin,outcomein,H0)
if ischar(H0)
    switch H0
        case 'Average'
            H0=mean(outcomein,'all','omitmissing');
        case 'Zero'
            H0=0;
    end
end
outcomein=repmat(outcomein',size(valsin,1),1);

group1=valsin .* outcomein;
valsout=nan(size(valsin,1),1);
psout=nan(size(valsin,1),1);
[~,psout,~,stats]=ttest(group1',H0);
psout=psout';
valsout = stats.tstat';
end