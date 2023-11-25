function [valsout,psout]=ea_explorer_stats_spearman(valsin,outcomein)
    [valsout,psout]=corr(valsin',outcomein,'Rows','pairwise','Type','Spearman');
end