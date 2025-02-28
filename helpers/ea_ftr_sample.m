function [fibers, idx] = ea_ftr_sample(ftrFile, sel, type, outputName)
% Sample fibers in the input FTR file
% 
% Arguments:
%     ftrFile: input FTR file
%     sel: number/interval/index of selected fibers
%     type: type of 'sel', can be 'number', 'interval' and 'index'
%     outputName: output file name
%
% Example:
%     ea_ftr_sample('FTR.mat', 20000, 'number');    Sample 20000 fibers
%     ea_ftr_sample('FTR.mat', 2, 'interval');      Sample with interval 2
%     ea_ftr_sample('FTR.mat', 1:2:20000, 'index'); Sample the specified fibers

if nargin < 2
    error('Please specify the sample indices!');
end

if isscalar(sel) % sel is a single number
    if nargin < 3 % Treat as sample number by default if type not specified
        warning('Number type not specified! Used as sample number.')
        type = 'number';
    end
else % sel is a serial of numbers
    type = 'index';
end

% load FTR file
ftr = load(ftrFile);

switch lower(type)
    case 'number'
        if sel < numel(ftr.idx)
            index = round(linspace(1, numel(ftr.idx), sel));
        else
            warning('Sample number larger than fiber number (%d)! Select all fibers.', numel(ftr.idx));
            index = 1:numel(ftr.idx);
        end
    case 'interval'
        index = 1:sel:numel(ftr.idx);
    case 'index'
        index = sel;
end

% Calculate selected fiber indices
endIndex = cumsum(ftr.idx);
startIndex = [1;endIndex+1];
startIndex = startIndex(1:end-1);
fibIndex = cell2mat(arrayfun(@(x,y) x:y, startIndex(index), endIndex(index), 'Uni', 0)');

% Sample fibers and idx
fibers = ftr.fibers(fibIndex,:);
idx = ftr.idx(index);

% Reset fibers(:, 4)
fibers(:, 4) = repelem(1:numel(idx), idx);

ftr.fibers = fibers;
ftr.idx = idx;

if exist('outputName', 'var')
    outputDir = fileparts(ftrFile);
    save(fullfile(outputDir, outputName), '-struct', 'ftr', '-v7.3');
end
