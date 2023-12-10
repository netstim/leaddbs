function [valsout,psout]=ea_explorer_stats_nmap(valsin)
valsout=sum(valsin,2,'omitmissing');
psout=zeros(size(valsout));
end