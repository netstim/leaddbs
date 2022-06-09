function [fibsvalBin, fibsvalSum, fibsvalMean, fibsvalPeak, fibsval5Peak, fibcell] = ea_discfibers_calcvals_cleartune(vatlist, fibcell, thresh)
% Calculate fiber connection values based on the VATs and the connectome

prefs = ea_prefs;
if ~exist('thresh','var')
    thresh = prefs.machine.vatsettings.horn_ethresh*1000;
end
if thresh < 150
    answer = questdlg("Threshold values less than 150 may cause errors in processing, setting default threshold value is 200");
    if strcmp(answer,'Yes')
        thresh=200;
    end
end
[numPatient, numSide] = size(vatlist);

fibsvalBin = cell(1, numSide);
fibsvalSum = cell(1, numSide);
fibsvalMean = cell(1, numSide);
fibsvalPeak = cell(1, numSide);
fibsval5Peak = cell(1, numSide);

for side = 1:numSide
    fibsvalBin{side} = zeros(length(fibcell{side}), numPatient);
    fibsvalSum{side} = zeros(length(fibcell{side}), numPatient);
    fibsvalMean{side} = zeros(length(fibcell{side}), numPatient);
    fibsvalPeak{side} = zeros(length(fibcell{side}), numPatient);
    fibsval5Peak{side} = zeros(length(fibcell{side}), numPatient);
    fibers=ea_fibcell2fibmat(fibcell{side});
    disp(['Calculate for side ', num2str(side), ':']);
    for pt = 1:numPatient
        disp(['VAT ', num2str(pt, ['%0',num2str(numel(num2str(numPatient))),'d']), '/', num2str(numPatient), '...']);
        if isstruct(vatlist) % direct nifti structs supplied
            vat=vatlist(pt,side);
        elseif iscell(vatlist) % filenames
            vat = ea_load_nii(vatlist{pt,side});
        end
        % Threshold the vat efield
        vatInd = find(abs(vat.img(:))>thresh);

        % Trim connectome fibers
        [xvox, yvox, zvox] = ind2sub(size(vat.img), vatInd);
        vatmm = ea_vox2mm([xvox, yvox, zvox], vat.mat);

        
        filter = all(fibers(:,1:3)>=min(vatmm),2) & all(fibers(:,1:3)<=max(vatmm), 2);

        % Skip further calculation in case VAT is totally not connected
        if ~any(filter)
            continue;
        end

        trimmedFiber = fibers(filter,:);

        % Map mm connectome fibers into VAT voxel space
        [trimmedFiberInd, ~, trimmedFiberID] = unique(trimmedFiber(:,4), 'stable');
        fibVoxInd = splitapply(@(fib) {ea_mm2uniqueVoxInd(fib, vat)}, trimmedFiber(:,1:3), trimmedFiberID);

        % Remove outliers
        fibVoxInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];
        trimmedFiberInd(cellfun(@(x) any(isnan(x)), fibVoxInd)) = [];

        % Find connected fibers
        connected = cellfun(@(fib) any(ismember(fib, vatInd)), fibVoxInd);

        % Generate binary fibsval for the T-test method
        fibsvalBin{side}(trimmedFiberInd(connected), pt)=1;

        % Check intersection between vat and the connected fibers
        vals = cellfun(@(fib) vat.img(intersect(fib, vatInd)), fibVoxInd(connected), 'Uni', 0);

        % Generate fibsval for the Spearman's correlation method
        fibsvalSum{side}(trimmedFiberInd(connected), pt) = cellfun(@sum, vals);
        fibsvalMean{side}(trimmedFiberInd(connected), pt) = cellfun(@mean, vals);
        fibsvalPeak{side}(trimmedFiberInd(connected), pt) = cellfun(@max, vals);
        fibsval5Peak{side}(trimmedFiberInd(connected), pt) = cellfun(@(x) mean(maxk(x,ceil(0.05*numel(x)))), vals);
    end

    % Keep all fibers in this case, convert to sparse matrix
    fibsvalBin{side} = sparse(fibsvalBin{side});
    fibsvalSum{side} = sparse(fibsvalSum{side});
    fibsvalMean{side} = sparse(fibsvalMean{side});
    fibsvalPeak{side} = sparse(fibsvalPeak{side});
    fibsval5Peak{side} = sparse(fibsval5Peak{side});
end


function fibers=ea_fibcell2fibmat(fibers)
[idx,~]=cellfun(@size,fibers);
fibers=cell2mat(fibers);
idxv=zeros(size(fibers,1),1);
lid=1; cnt=1;
for id=idx'

    idxv(lid:lid+id-1)=cnt;
    lid=lid+id;
    cnt=cnt+1;
end
fibers=[fibers,idxv];
