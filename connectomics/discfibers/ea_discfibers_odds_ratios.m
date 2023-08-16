function [OR,CI95_up,CI95_low,p] = ea_discfibers_odds_ratios(activ_prob, outcome)

% See the method's description in Jergas et al 2023 https://doi.org/10.1101/2023.04.26.23289100


%disp('Calculating fiber-wise Odds Ratios for EF-sigmoids ...')

allAntiValues = 1 - activ_prob; %%non-connected voxvals = 1-voxval
for f = 1:size(activ_prob,2)
    %     OR_table(f,1) = sum(convals(f,find(dys == 1 & conmatrix(f,:) == 1))); %dys_stim
    %     OR_table(f,2) = sum(inverse_convals(f,find(dys == 1 & conmatrix(f,:) == 0))); %dys_no_stim
    %     OR_table(f,3) = sum(convals(f,find(dys == 0 & conmatrix(f,:) == 1))); %no_dys_stim
    %     OR_table(f,4) = sum(inverse_convals(f,find(dys == 0 & conmatrix(f,:) == 0))); %no_dys_no_stim
    OR_table(f,1) = sum(activ_prob(outcome,f)); %dys_stim
    OR_table(f,2) = sum(allAntiValues(outcome,f)); %dys_no_stim
    OR_table(f,3) = sum(activ_prob(~outcome,f)); %no_dys_stim
    OR_table(f,4) = sum(allAntiValues(~outcome,f)); %no_dys_no_stim
    if any(OR_table(f,:) == 0)
        OR_table(f,:) = OR_table(f,:) + 0.5; %% Haldane-Anscombe correction s. Zero-cell corrections in random-effects meta-analyses; Frank Weber*, Guido Knapp†, Katja Ickstadt†, Günther Kundt*, Änne Glass
    end
end


for f = 1:height(OR_table)
    OR(f) = (OR_table(f,1) * OR_table(f,4)) / (OR_table(f,2) * OR_table(f,3));
    CI95_up(f) =  exp(log(OR(f))+1.96*sqrt(1/OR_table(f,1)+1/OR_table(f,2)+1/OR_table(f,3)+1/OR_table(f,4))); %%https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2938757/
    CI95_low(f) = exp(log(OR(f))-1.96*sqrt(1/OR_table(f,1)+1/OR_table(f,2)+1/OR_table(f,3)+1/OR_table(f,4)));
    Est = log(OR(f));
    SE = (log(CI95_up(f))-log(CI95_low(f)))/(2*1.96);
    z = Est/SE;
    if z<0
        z = -1*z;
    end
    p(f) = exp((-0.717*z)-(0.416*z^2)); %%https://www.bmj.com/content/343/bmj.d2304
end
