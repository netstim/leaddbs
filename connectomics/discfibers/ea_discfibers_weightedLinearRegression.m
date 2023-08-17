function [vals,p_vals] = ea_discfibers_weightedLinearRegression(activ_prob, outcome, expected_outcome_metric, output_stat)

% by Till Dembek, publication TBA
% Warning: this metric takes time to be calculated
        
expected_outcome = zeros(size(activ_prob));
if strcmp(expected_outcome_metric,'mean')
    expected_outcome(:) = ea_nanmean(outcome);
end

outcome_over_fibers = repmat(outcome,[1,size(activ_prob,2)]);

vals = zeros(size(activ_prob,2),1);
p_vals = zeros(size(activ_prob,2),1);

mysyntax = 'outcome ~ 1 + condition';
for l=1:size(activ_prob,2)

    %% wilcoxon signed rank test for every Voxel in VoxelSum versus VoxelMean
    mytable = table;
    mytable.condition = vertcat(zeros(size(expected_outcome(:,l))),ones(size(outcome_over_fibers(:,l))));
    mytable.outcome = vertcat(expected_outcome(:,l),outcome_over_fibers(:,l));
    mytable.weight = vertcat(activ_prob(:,l),activ_prob(:,l));
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
