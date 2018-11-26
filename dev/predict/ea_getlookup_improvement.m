function Ihat=ea_getlookup_improvement(ptvta,othervta,otherI,exponent)
if ~exist('exponent','var')
    exponent=1;
end
weights=sum(othervta.*repmat(ptvta,size(othervta,1),1),2);
weights=(weights).^exponent;
weights=weights./sum(weights);
Ihat=sum(otherI.*weights);











