function ea_corrfib2trk(reformatted_corrfibs, sel, outputName)
% Convert correlative fibertracts to trk file

if nargin < 2 || isempty(sel)
    sel = 'both';
end

load(reformatted_corrfibs);

% Select fibers
if ischar(sel)
    switch lower(sel)
        case 'both'
            selInd = 1:length(fibcell);
            suffix = '_all';
        case {'positive', 'pos'}
            selInd = vals > 0;
            suffix = '_pos';
        case {'negative', 'neg'}
            selInd = vals < 0;
            suffix = '_neg';
    end
else % 'sel' is a threshold
    selInd = vals > sel;
    suffix = ['_th', num2str(sel)];
end

if nargin < 3
    outputName = ['corrFTR', suffix];
else
    outputName = strrep(outputName, '.mat', '');
end

fibcell = fibcell(selInd);

% Construct fibers for trk conversion
idx = cellfun(@length, fibcell);
fibers = zeros(sum(idx), 4);
fibers(:, 1:3) = cell2mat(fibcell);
fibInd = repelem((1:length(idx))', idx);
fibers(:, 4) = fibInd;

% Meta information
ea_fibformat = '1.1';
fourindex = 1;
voxmm = 'mm';

% Save FTR for trk conversion
if isempty(fileparts(reformatted_corrfibs))
    outputDir = '.';
else
    outputDir = fileparts(reformatted_corrfibs);
end
save(fullfile(outputDir, outputName), ...
     'ea_fibformat', 'fibers', 'fourindex', 'idx', 'voxmm', '-v7.3');

% FTR to TRK conversion
ea_ftr2trk(fullfile(outputDir, outputName));
