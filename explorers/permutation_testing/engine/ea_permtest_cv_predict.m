function [out, firstModel] = ea_permtest_cv_predict(cvb, O, iters)
% k-fold cross-validated predictions of the old fiber filtering explorer,
% for several versions of the scores at once (e.g. the real scores and a
% batch of shuffles). Repeats, step by step, ea_disctract.kfoldcv ->
% crossval -> ea_compute_fibscore_model -> ea_discfibers_calcstats:
%
% for every k-fold iteration and fold, on the fold's training patients:
%   1. coverage filter (fibers reached by too few E-fields are dropped)
%   2. fiber test of every model (subscore) -- the fast path computes it
%      for all versions with one matrix product (ea_permtest_fastcorr*.m)
%   3. optional significance masking over all models and sides (ea_corrsignan)
%   4. display threshold per model and side (copy of ea_fibValThresh)
%   5. optional model normalization (as the explorer: vals{s} for
%      s = 1:length(vals), i.e. in subscore mode the right side of each
%      model only)
%   6. score the test patients: spatial correlation / mean / sum / peak of
%      the model with their (unmirrored) fiber values; patients without any
%      fiber value on the model get NaN
% then per iteration
%   7. optional fit to scores (fitglm per fold and model on training scores)
%   8. subscore mode: weighted combination of the models (subscore weights)
% and finally the mean over iterations.
% A model without fibers in a fold predicts NaN for that fold's test patients.
%
% cvb   - from ea_permtest_cv_adapter_disctract.m.
% O     - Npatients x (Nvoters + 2) x B scores per version: columns
%         1:Nvoters = models' scores, Nvoters+1 = combined target (NaN if
%         none), Nvoters+2.. = weight variables (one per model, NaN if fixed).
% iters - k-fold iterations to use (default all; e.g. 1 for a timing probe).
%
% out - struct:
%   .voterPred   Nsel x Nvoters x B, prediction of every model (after fit to
%                scores if on), mean over iterations
%   .combPred    Nsel x B, combined prediction (subscore mode) or the single
%                model's prediction
%   .nModelFibers  numel(iters) x kfold x Nvoters x Nsides x B, fibers in each fold's model
%   .nAnalysed     numel(iters) x kfold x Nsides, fibers after the coverage filter
%   .fast          true when every fiber test used the fast path
% firstModel - optional: the model of the first iteration, first fold and
%              first version (.vals, .usedidx as rows of cvb.fib; for checks)

cv = cvb.cv;
if ~exist('iters', 'var') || isempty(iters)
    iters = 1:cv.kIter;
end
N = cvb.Npatients;
sel = cvb.patientselection(:)';
Nsel = numel(sel);
V = numel(cvb.voterLabels);
B = size(O, 3);
Nsides = cvb.Nsides;
connthr = cvb.statsettings.connthreshold / 100;
efthr = cvb.statsettings.efieldthreshold;
nrows = cellfun(@(x) size(x, 1), cvb.fib.model);
multcomp = struct('multcompstrategy', cv.multcompstrategy, 'alphalevel', cv.alphalevel);

voterPredIt = nan(Nsel, V, B, numel(iters));
combPredIt = nan(Nsel, B, numel(iters));
nModel = zeros(numel(iters), cv.kfold, V, Nsides, B, 'uint32');
nAnal = zeros(numel(iters), cv.kfold, Nsides, 'uint32');
allFast = true;

for ii = 1:numel(iters)
    it = iters(ii);
    Ihat = nan(Nsel, V, B);
    if cv.doactualprediction
        Itrain = nan(cv.kfold, Nsel, V, B);
    end
    for c = 1:cv.kfold
        training = cv.folds(:, it)' ~= c;
        test = ~training;
        trsel = sel(training);
        gpatsel = trsel;
        if cvb.mirrorsides
            gpatsel = [trsel, trsel + N];
        end
        if cv.doactualprediction
            cols = 1:Nsel;
        else
            cols = find(test);
        end
        testInCols = test(cols);

        % 1. coverage filter on the training patients
        nonempty = cell(1, Nsides);
        valsin = cell(1, Nsides);
        for side = 1:Nsides
            G = cvb.fib.model{side}(:, gpatsel);
            Nmap = sum(G > efthr, 2);
            keep = Nmap >= connthr * numel(gpatsel);
            nonempty{side} = find(keep & sum(G, 2, 'omitnan') > 0);
            valsin{side} = G(nonempty{side}, :);
            nAnal(ii, c, side) = numel(nonempty{side});
        end

        % 2. fiber test of every model, all versions at once
        R = cell(V, Nsides);
        P = cell(V, Nsides);
        for v = 1:V
            I = reshape(O(:, v, :), N, B);
            if cvb.mirrorsides
                I = [I; I]; %#ok<AGROW>
            end
            Y = I(gpatsel, :);
            for side = 1:Nsides
                if isempty(nonempty{side})
                    R{v, side} = zeros(0, B);
                    P{v, side} = zeros(0, B);
                    continue
                end
                prep = ea_permtest_fastcorr_prepare(cvb.statfile, valsin{side}, Y(:, 1));
                if prep.ok
                    [R{v, side}, P{v, side}] = ea_permtest_fastcorr(prep, Y);
                else
                    allFast = false;
                    R{v, side} = nan(numel(nonempty{side}), B);
                    P{v, side} = nan(numel(nonempty{side}), B);
                    for b = 1:B
                        [r, p] = feval(cvb.statfile, valsin{side}, Y(:, b), cvb.H0);
                        R{v, side}(:, b) = r(:);
                        P{v, side}(:, b) = p(:);
                    end
                end
            end
        end

        % 3.-6. per version
        for b = 1:B
            vals = cell(V, Nsides);
            pvals = cell(V, Nsides);
            for v = 1:V
                for side = 1:Nsides
                    vals{v, side} = nan(nrows(side), 1);
                    vals{v, side}(nonempty{side}) = R{v, side}(:, b);
                    pvals{v, side} = nan(nrows(side), 1);
                    pvals{v, side}(nonempty{side}) = P{v, side}(:, b);
                end
            end
            if cv.showsignificantonly
                vals = ea_corrsignan(vals, pvals, multcomp);
            end
            usedidx = cell(V, Nsides);
            for v = 1:V
                for side = 1:Nsides
                    usedidx{v, side} = find(isfinite(vals{v, side}));
                    vals{v, side} = vals{v, side}(usedidx{v, side});
                end
            end
            [vals, usedidx] = local_threshold(vals, usedidx, cv);
            switch cv.modelNormalization
                case 'z-score'
                    for s = 1:length(vals)
                        vals{s} = ea_nanzscore(vals{s});
                    end
                case 'van Albada 2007'
                    for s = 1:length(vals)
                        vals{s} = ea_normal(vals{s});
                    end
            end
            if nargout > 1 && ii == 1 && c == 1 && b == 1
                firstModel = struct('vals', {vals}, 'usedidx', {usedidx});
            end
            for v = 1:V
                for side = 1:Nsides
                    nModel(ii, c, v, side, b) = numel(usedidx{v, side});
                end
                if all(cellfun(@isempty, vals(v, :)))
                    continue
                end
                vflat = vertcat(vals{v, :});
                F = cell(Nsides, 1);
                for side = 1:Nsides
                    F{side} = cvb.fib.pred{side}(usedidx{v, side}, sel(cols));
                end
                F = vertcat(F{:});
                ih = local_score(vflat, F, cv);
                ih(testInCols & ~sum(F, 1)) = NaN; % no fiber value on the model
                Ihat(test, v, b) = ih(testInCols);
                if cv.doactualprediction
                    Itrain(c, training, v, b) = ih(~testInCols);
                end
            end
        end
    end

    % 7. fit to scores (per fold, on the fold's training patients)
    if cv.doactualprediction
        fitted = nan(Nsel, V, B);
        for b = 1:B
            for v = 1:V
                useI = O(sel, v, b);
                for c = 1:cv.kfold
                    training = cv.folds(:, it) ~= c;
                    test = ~training;
                    mdl = fitglm(reshape(Itrain(c, training, v, b), [], 1), useI(training), cv.predictionmodel);
                    fitted(test, v, b) = predict(mdl, Ihat(test, v, b));
                end
            end
        end
        Ihat = fitted;
    end
    voterPredIt(:, :, :, ii) = Ihat;

    % 8. combination of the models
    if cv.subscoreMode
        for b = 1:B
            W = zeros(Nsel, V);
            for v = 1:V
                if ~isnan(cv.weights(v))
                    W(:, v) = cv.weights(v);
                else
                    w = O(sel, V + 1 + v, b);
                    W(:, v) = ea_minmax(w) ./ max(w);
                end
            end
            for i = 1:Nsel
                W(i, :) = W(i, :) ./ ea_nansum(W(i, :));
            end
            combPredIt(:, b, ii) = ea_nansum(Ihat(:, :, b) .* W, 2);
        end
    else
        combPredIt(:, :, ii) = reshape(Ihat(:, 1, :), Nsel, B);
    end
end

out = struct;
out.voterPred = mean(voterPredIt, 4, 'omitnan');
out.combPred = mean(combPredIt, 3, 'omitnan');
out.nModelFibers = nModel;
out.nAnalysed = nAnal;
out.fast = allFast;


function [vals, usedidx] = local_threshold(vals, usedidx, cv)
% Display threshold per model and side, as in ea_discfibers_calcstats.m
unthresholded = vals;
[V, Nsides] = size(vals);
for v = 1:V
    for side = 1:Nsides
        if strcmp(cv.threshstrategy, 'Fixed Amount') % each side separately
            posvals = sort(vals{v, side}(vals{v, side} > 0), 'descend');
            negvals = sort(vals{v, side}(vals{v, side} < 0), 'ascend');
        else % the same threshold for both sides
            allvals = vertcat(unthresholded{v, :});
            posvals = sort(allvals(allvals > 0), 'descend');
            negvals = sort(allvals(allvals < 0), 'ascend');
        end
        if ~cv.vis.pos(v) || ~cv.vis.posamount(v, side) || isempty(posvals)
            posthresh = inf;
        else
            posthresh = local_fibValThresh(cv.threshstrategy, posvals, cv.vis.posamount(v, side));
        end
        if ~cv.vis.neg(v) || ~cv.vis.negamount(v, side) || isempty(negvals)
            negthresh = -inf;
        else
            negthresh = local_fibValThresh(cv.threshstrategy, negvals, cv.vis.negamount(v, side));
        end
        remove = vals{v, side} < posthresh & vals{v, side} > negthresh;
        vals{v, side}(remove) = [];
        usedidx{v, side}(remove) = [];
    end
end


function fibValThreshold = local_fibValThresh(threshstrategy, vals, threshold)
% Copy of ea_fibValThresh (local function in ea_discfibers_calcstats.m)
switch threshstrategy
    case 'Percentage Relative to Peak'
        range = vals(1) - vals(end);
        fibValThreshold = vals(1) - threshold/100 * range;
        if range == 0
            if vals(1) > 0
                fibValThreshold = fibValThreshold - eps*10;
            else
                fibValThreshold = fibValThreshold + eps*10;
            end
        end
    case 'Percentage Relative to Amount'
        index = round((threshold/100)*length(vals));
        if index <=0
            fibValThreshold = vals(1);
        else
            fibValThreshold = vals(index);
        end
    case 'Fixed Amount'
        if length(vals)>round(threshold)
            fibValThreshold=vals(round(threshold));
        else
            fibValThreshold=vals(end);
        end
    case 'Histogram (CDF)'
        if vals(1) > 0
            [fx, x] = ecdf(vals);
            fibValThreshold = x(find(fx>=(1-threshold), 1));
        else
            [fx, x] = ecdf(-vals);
            fibValThreshold = -x(find(fx>=(1-threshold), 1));
        end
    case 'Fixed Fiber Value'
        fibValThreshold = threshold;
end


function ih = local_score(vflat, F, cv)
% Score patients (columns of F) with the model vflat, as
% ea_compute_fibscore_model.m does for E-fields (both hemispheres together)
switch cv.basepredictionon
    case 'profile of scores: spearman'
        ih = corr(vflat, F, 'rows', 'pairwise', 'type', 'spearman');
    case 'profile of scores: pearson'
        ih = corr(vflat, F, 'rows', 'pairwise', 'type', 'pearson');
    case 'profile of scores: bend'
        ih = ea_bendcorr(vflat, F);
    case 'mean of scores'
        ih = ea_nanmean(vflat .* F, 1);
    case 'sum of scores'
        ih = ea_nansum(vflat .* F, 1);
    case 'peak of scores'
        ih = ea_discfibers_getpeak(vflat .* F, cv.peakvis(1), cv.peakvis(2), 'peak');
    case 'peak 5% of scores'
        ih = ea_discfibers_getpeak(vflat .* F, cv.peakvis(1), cv.peakvis(2), 'peak5');
end
ih = reshape(ih, 1, []);
