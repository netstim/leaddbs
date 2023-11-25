function [valsout,psout]=ea_explorer_stats_pearson(valsin,outcomein)
[valsout,psout]=corr(valsin',outcomein,'Rows','pairwise','Type','Pearson');
end