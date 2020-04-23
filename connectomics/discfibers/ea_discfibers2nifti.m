function ea_discfibers2nifti(discfiber, threshold, fiberset, outputName, reference)
% Convert discriminative fibers to nifti
%
%     discfiber: discriminative fibers including fiber cell and T-scores
%     threshold: threshold to filter the fibers
%     fiberset: fiber set to be converted, can be 'postive', 'negative' or 'both'
%     outputName: output file name
%     reference: reference nifti which defines the space

load(discfiber, 'fibcell')
load(discfiber, 'vals')

if ~exist('threshold', 'var')
    threshold = 5/100;
elseif isempty(threshold)
    threshold = 0.05;
elseif threshold>1  % Use percentage
    threshold = threshold/100;
end

if ~exist('fiberset', 'var')
    fiberset = 'pos';
end

if ~exist('outputName', 'var')
    outputName = regexprep(discfiber, '\.mat$', '.nii');
end

if ~exist('reference', 'var')
    reference = [ea_space, 't1.nii'];
end

vals(isnan(vals))=0;
posits = vals(vals>0);
posits = sort(posits,'descend');
negits = vals(vals<0);
negits = sort(negits,'ascend');
% Determine positive/negative threshold
if ismember(fiberset, {'pos', 'positive'})
    posthresh = posits(round(length(posits)*threshold));
    negthresh = min(vals)-eps;
    disp(['Fiber colors: Positive (T = ',num2str(posthresh),' ~ ',num2str(posits(1)), ')']);
elseif ismember(fiberset, {'neg', 'negative'})
    posthresh = max(vals)+eps;
    negthresh = negits(round(length(negits)*threshold));
    disp(['Fiber colors: Negative (T = ',num2str(negits(1)),' ~ ',num2str(negthresh), ')']);
elseif strcmp(fiberset, 'both')
    posthresh = posits(round(length(posits)*threshold));
    negthresh = negits(round(length(negits)*threshold));
    disp(['Fiber colors: Positive (T = ',num2str(posthresh),' ~ ',num2str(posits(1)), ')']);
    disp(['Fiber colors: Negative (T = ',num2str(negits(1)),' ~ ',num2str(negthresh), ')']);
end

% Remove vals and fibers outside the thresholding range
remove = logical(logical(vals<posthresh) .* logical(vals>negthresh));
vals(remove)=[];
fibcell(remove)=[];

% Convert fibers from mm to voxel
refnii = ea_load_nii(reference);
fibcell_vox = cellfun(@(x) round(ea_mm2vox(x, refnii.mat)), fibcell, 'Uni', 0);

% Remove negative coordinates
keep = cellfun(@(x) sum(x>0,2)==3, fibcell_vox, 'Uni', 0);
fibcell = cellfun(@(x, y) x(y,:), fibcell, keep, 'Uni', 0);
fibcell_vox = cellfun(@(x, y) x(y,:), fibcell_vox, keep, 'Uni', 0);

% Remove single point cell
keep = cellfun(@(x) size(x,1), fibcell) > 1;
fibcell = fibcell(keep);
fibcell_vox = fibcell_vox(keep);
vals = vals(keep);

% Calculate fiber indices
idx = cellfun(@(x) size(x,1), fibcell);
fiberval = repelem(vals, idx);

% Convert fibcell to fiber matrix
fibers_vox = cell2mat(fibcell_vox);

% Calculate values at voxels
refnii.img = zeros(size(refnii.img));
fibImgInd = sub2ind(size(refnii.img), fibers_vox(:,1), fibers_vox(:,2), fibers_vox(:,3));
[unifibImgInd, ~, ic] = unique(fibImgInd);
fibImgVal = accumarray(ic, fiberval);
refnii.img(unifibImgInd) = refnii.img(unifibImgInd) + fibImgVal;

% Invert the values when exporting negative fiberset
if ismember(fiberset, {'neg', 'negative'})
    refnii.img = refnii.img*-1;
end

% Save fiber nifti
refnii.fname = outputName;
ea_write_nii(refnii);
