function [valsout,psout]=ea_explorer_stats_signedranktest(valsin,outcomein,H0)
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
psout=nan(size(valsin,1),1);
for i=1:size(valsin,1)
    [psout(i),~,stats(i)]=signrank(group1(i,:),repmat(H0,size(group1(i,:))),'method','exact');
    testsign(i,1) =  sign(median(group1(i,:),"all",'omitmissing')-H0);
end
valsout = testsign.*-log10(psout);
end