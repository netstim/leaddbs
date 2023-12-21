function [vals,p_vals] = ea_discfibers_TwoSample_weightedLinearRegression(activ_prob, outcome, output_stat)

% by Till Dembek, publication TBA
% Warning: this metric takes time to be calculated
       
outcome_over_fibers = repmat(outcome,[1,size(activ_prob,2)]);

vals = zeros(size(activ_prob,2),1);
p_vals = zeros(size(activ_prob,2),1);

allAntiValues = 1 - activ_prob;

mysyntax = 'outcome ~ 1 + condition';
for l=1:size(activ_prob,2)

    mytable = table;
    mytable.condition = vertcat(zeros(size(activ_prob(:,l))),ones(size(allAntiValues(:,l))));
    mytable.outcome = vertcat(outcome_over_fibers(:,l),outcome_over_fibers(:,l));
    mytable.weight = vertcat(allAntiValues(:,l),activ_prob(:,l));
    mymdl = fitlm(mytable,mysyntax,'Weights',mytable.weight);

    if strcmp(output_stat,'t-value')
        vals(l) = mymdl.Coefficients.tStat(2);
    elseif strcmp(output_stat,'slop')
        vals(l) = mymdl.Coefficients.Estimate(2);
    end
    p_vals(l) = mymdl.Coefficients.pValue(2);


    clear mymdl mytable
end

end
