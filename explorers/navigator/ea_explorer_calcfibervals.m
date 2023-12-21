function [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell, connFiberInd, totalFibers] = ea_explorer_calcfibervals(vatlist, cfile, thresh)
% Calculate fiber connection values based on the VATs and the connectome

disp('Load Connectome...');
load(cfile, 'fibers', 'idx');

prefs = ea_prefs;
if ~exist('thresh','var')
    thresh = prefs.machine.vatsettings.horn_ethresh*1000;
end
[numPatient, numSide] = size(vatlist);

fibsvalBin = cell(1, numSide);
fibsvalSum = cell(1, numSide);
fibsvalMean = cell(1, numSide);
fibsvalPeak = cell(1, numSide);
fibsval5Peak = cell(1, numSide);

fibcell = cell(1, numSide);
connFiberInd = cell(1, numSide);

totalFibers = length(idx); % total number of fibers in the connectome to work with global indices


for side = 1:numSide
    %% crop fibers to template
    fibsvalBin{side} = nan(totalFibers, numPatient);
    fibsvalSum{side} = nan(totalFibers, numPatient);
    fibsvalMean{side} = nan(totalFibers, numPatient);
    fibsvalPeak{side} = nan(totalFibers, numPatient);
    fibsval5Peak{side} = nan(totalFibers, numPatient);

    % disp(['Calculate for side ', num2str(side), ':']);
    ea_dispercent(1/numPatient,['Calculate fiber Activation for side ', num2str(side)])
    for pt = 1:numPatient
        if isstruct(vatlist) % direct nifti structs supplied
            vat = vatlist(pt,side);
        elseif iscell(vatlist) % filenames
            if isfile(vatlist{pt,side})
                vat = ea_load_nii(vatlist{pt,side});
            else
                ea_cprintf('CmdWinWarnings', 'Skipping calculating connectivity: VTA doesn''t exist!\n');
                continue;
            end
        end
        % Threshold the vat efield
        vat.img(vat.img<thresh) = nan;
        [fibidx,peakvals,meanvals,sumvals,peak5vals,binvals] = ea_explorer_efoverlapfiber(vat,vat.img,fibers);
        % Generate fibsval for the Spearman's correlation method
        fibsvalBin{side}(fibidx,pt) = binvals;
        fibsvalSum{side}(fibidx,pt) = sumvals;
        fibsvalMean{side}(fibidx,pt) = meanvals;
        fibsvalPeak{side}(fibidx,pt) = peakvals;
        fibsval5Peak{side}(fibidx,pt) = peak5vals;
        clear fibidx peakvals meanvals sumvals binvals peak5vals
        ea_dispercent(pt/numPatient)
    end
    ea_dispercent(pt/numPatient,'end')
    disp('Completed.')
    % Remove values for not connected fibers, convert to sparse matrix
    fibIsConnected = any(fibsvalBin{side}, 2);
    fibsvalBin{side} = sparse(fibsvalBin{side}(fibIsConnected, :));
    fibsvalSum{side} = sparse(fibsvalSum{side}(fibIsConnected, :));
    fibsvalMean{side} = sparse(fibsvalMean{side}(fibIsConnected, :));
    fibsvalPeak{side} = sparse(fibsvalPeak{side}(fibIsConnected, :));
    % fibsval5Peak{side} = sparse(fibsval5Peak{side}(fibIsConnected, :));

    % Extract connected fiber cell
    connFiberInd{side} = find(fibIsConnected);
    connFiber = fibers(ismember(fibers(:,4), connFiberInd{side}), 1:3);
    fibcell{side} = mat2cell(connFiber, idx(connFiberInd{side}));
end



