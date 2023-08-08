function ea_discfibers2trk(disfiber, sel, outputName, LPS)
% Convert discriminitive fibertracts to trk file
%
% If the trk is going to be visualized in Surf-Ice, LPS should be set to 1
% to fix the orientation.

if nargin < 2 || isempty(sel)
    sel = 'both';
end

load(disfiber, 'fibcell')
load(disfiber, 'vals')

if size(fibcell,2) == 2
    fibcell = vertcat(fibcell{:});
end

if iscell(vals) && size(vals,2) == 2
    vals = vertcat(vals{:});
end

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
    outputName = strrep(outputName, '.trk', '');
end

fibcell = fibcell(selInd);

% Construct fibers for trk conversion
idx = cellfun(@(x) size(x,1), fibcell);
fibers = zeros(sum(idx), 4);
fibers(:, 1:3) = cell2mat(fibcell);
fibInd = repelem((1:length(idx))', idx);
fibers(:, 4) = fibInd;
vals = vals(selInd);

% Meta information
ea_fibformat = '1.1';
fourindex = 1;
voxmm = 'mm';

% Save FTR for trk conversion
if isempty(fileparts(disfiber))
    outputDir = '.';
else
    outputDir = fileparts(disfiber);
end
save(fullfile(outputDir, [outputName, '.mat']), ...
     'ea_fibformat', 'fibers', 'fourindex', 'idx', 'voxmm', 'vals', '-v7.3');

% FTR to TRK conversion
ea_ftr2trk(fullfile(outputDir, [outputName, '.mat']), [], LPS);
