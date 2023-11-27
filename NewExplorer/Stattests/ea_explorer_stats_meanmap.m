function [valsout,psout]=ea_explorer_stats_meanmap(valsin,outcomein)
outcomein=repmat(outcomein',size(valsin,1),1);
valsin=~isnan(valsin); % valsin already only includes values above the threshold;
valsout=sum((outcomein.*valsin),2)./sum(valsin,2);
psout=zeros(size(valsout));
end