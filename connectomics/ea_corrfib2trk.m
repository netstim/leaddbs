function ea_corrfib2trk(reformatted_corrfibs, sel)
% Convert correlative fibertracts to trk file

if nargin < 2
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
save(fullfile(fileparts(reformatted_corrfibs), ['corrFTR', suffix, '.mat']), ...
     'ea_fibformat', 'fibers', 'fourindex', 'idx', 'voxmm', '-v7.3');

% FTR to TRK conversion
ea_ftr2trk(['corrFTR', suffix], fileparts(reformatted_corrfibs));
