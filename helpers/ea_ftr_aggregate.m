function [fibers, idx] = ea_ftr_aggregate(ftrFiles, outputFile, sel, type)
% Aggregate fibers in the input FTR files
% 
% Arguments:
%     ftrFiles: input FTR files
%     outputFile: output file
%     sel: number/interval/index of selected fibers
%     type: type of 'sel', can be 'number', 'interval' and 'index'
%
% Example:
%     ea_ftr_aggregate(ftrFiles, 'FTR.mat');                     Aggregate fibers in ftrFiles
%     ea_ftr_aggregate(ftrFiles, 'FTR.mat', 20000, 'number');    Sample 20000 fibers and aggregate
%     ea_ftr_aggregate(ftrFiles, 'FTR.mat', 2, 'interval');      Sample with interval 2 and aggregate
%     ea_ftr_aggregate(ftrFiles, 'FTR.mat', 1:2:20000, 'index'); Sample the specified fibers and aggregate

fibers = [];
idx = [];
offset = 0;

for i=1:length(ftrFiles)
    if nargin < 3 % Simple aggregate
        ftr = load(ftrFiles{i});
    else % Sample and aggregate
        if numel(sel) == 1 % sel is a single number
            if nargin < 4 % Treat as sample number by default if type not specified
                warning('Number type not specified! Used as sample number.')
                type = 'number';
            end
        else % sel is a serial of numbers
            type = 'index';
        end
        [ftr.fibers, ftr.idx] = ea_ftr_sample(ftrFiles{i}, sel, type);
    end
    
    ftr.fibers(:,4) = ftr.fibers(:,4) + offset;
    fibers = [fibers;ftr.fibers];
    idx = [idx;ftr.idx]; 
    offset = offset + numel(idx);
end

% Use other meta info from the first FTR file
ftr = load(ftrFiles{1});
ftr.fibers = fibers;
ftr.idx = idx;
save(outputFile, '-struct', 'ftr', '-v7.3');
