function [valsout,psout]=ea_explorer_stats_proportiontest(valsin,outcomein)
valsout=nan(size(valsin,1),1);
psout=nan(size(valsin,1),1);
% outcomein=round(outcomein./max(outcomein));
if ~all(outcomein==1|outcomein==0|isnan(outcomein))
    warning('off','backtrace')
    warning('Values other than 0 or 1 detected. Outcome is not binary. Please choose different test!')
    warning('on','backtrace')
    return
else
    numtrue = sum(outcomein==1);
    numfalse = sum(outcomein==0);
    
    outcomein=repmat(outcomein',size(valsin,1),1);
    valsin=valsin .* outcomein;
    sumtrue = sum(valsin==1,2);
    sumfalse = sum(valsin==0,2);

    for i=1:size(valsin,1)
        [~,psout(i), valsout(i)]  = ea_prop_test([sumtrue(i),sumfalse(i)],[numtrue,numfalse],1);
    end
end
end