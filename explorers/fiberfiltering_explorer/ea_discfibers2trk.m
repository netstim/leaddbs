function ea_discfibers2trk(discfiber, fiberset, outputName, LPS)
% Convert discriminitive fibertracts to trk file
%
% If the trk is going to be visualized in Surf-Ice, LPS should be set to 1
% to fix the orientation.

arguments
    discfiber {mustBeTextScalar}
    fiberset {mustBeMember(fiberset, {'both', 'pos', 'positive', 'neg', 'negative'})} = 'both'
    outputName {mustBeTextScalar} = ''
    LPS {mustBeNumericOrLogical} = false
end

load(discfiber, 'fibcell')
load(discfiber, 'vals')

if size(fibcell,2) == 2
    fibcell = vertcat(fibcell{:});
elseif size(fibcell,2) == 1
    fibcell = fibcell{:};
end

if iscell(vals) && size(vals,2) == 2
    vals = vertcat(vals{:});
elseif iscell(vals) && size(vals,2) == 1
    vals = vals{:};
end

% Select fibers
if ischar(fiberset)
    switch lower(fiberset)
        case 'both'
            selectedInd = 1:length(fibcell);
            suffix = '_posneg';
        case {'positive', 'pos'}
            selectedInd = vals > 0;
            suffix = '_pos';
        case {'negative', 'neg'}
            selectedInd = vals < 0;
            suffix = '_neg';
    end
else % 'sel' is a threshold
    selectedInd = vals > fiberset;
    suffix = ['_th', num2str(fiberset)];
end

if isempty(outputName)
    outputName = regexprep(discfiber, '\.mat$', suffix);
else
    outputName = strrep(outputName, '.trk', '');
end

fibcell = fibcell(selectedInd);

% Construct fibers for trk conversion
idx = cellfun(@(x) size(x,1), fibcell);
fibers = zeros(sum(idx), 4);
fibers(:, 1:3) = cell2mat(fibcell);
fibInd = repelem((1:length(idx))', idx);
fibers(:, 4) = fibInd;
vals = vals(selectedInd);

% Meta information
ea_fibformat = '1.1';
fourindex = 1;
voxmm = 'mm';

% Save FTR for trk conversion
if isempty(fileparts(discfiber))
    outputDir = '.';
else
    outputDir = fileparts(discfiber);
end
save(fullfile(outputDir, [outputName, '.mat']), ...
     'ea_fibformat', 'fibers', 'fourindex', 'idx', 'voxmm', 'vals', '-v7.3');

% FTR to TRK conversion
ea_ftr2trk(fullfile(outputDir, [outputName, '.mat']), [], LPS);
