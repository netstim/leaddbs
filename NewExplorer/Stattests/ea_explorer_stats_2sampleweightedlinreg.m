function [valsout,psout]=ea_explorer_stats_2sampleweightedlinreg(valsin,outcomein)

outcomein=repmat(outcomein',size(valsin,1),1);

group1=outcomein;
group1(isnan(valsin))=NaN;
group2=outcomein;
group2(~isnan(valsin))=NaN;
valsout=nan(size(valsin,1),1);
psout=nan(size(valsin,1),1);

mysyntax = 'outcome ~ 1+condition';
ea_dispercent(1/size(valsin,1),'Calculating Weighted regression')
for i = 1:size(valsin,1)
    ea_dispercent(i/size(valsin,1))
    mytable = table;
    mytable.outcome= vertcat(group2(i,:)',group1(i,:)');
    mytable.weight = vertcat(valsin(i,:)',1-valsin(i,:)');
    mytable.condition = vertcat(zeros(size(valsin(i,:)')),ones(size(valsin(i,:)')));
    mymdl = fitlm(mytable,mysyntax,'Weights',mytable.weight);
    %% storing statitic values in maps
    psout(i) = mymdl.Coefficients.pValue(2);
    valsout(i) = mymdl.Coefficients.tStat(2);
end
ea_dispercent(i/size(valsin,1),'end')
end