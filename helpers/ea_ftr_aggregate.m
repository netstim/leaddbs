function [fibers, idx] = ea_ftr_aggregate(ftrFiles, outputFile, sel, type, filtermask)
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

if ~exist('filtermask','var')
    filtermask='';
end

fibers = [];
idx = [];

for i=1:length(ftrFiles)
    fprintf('Aggregating fibers from %d/%d FTR files...\n', i, length(ftrFiles));
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

    fibers = [fibers;ftr.fibers];
    idx = [idx;ftr.idx]; 
end

% Reset fibers(:, 4)
fibers(:, 4) = repelem(1:numel(idx), idx);

% filter
if ~isempty(filtermask)
    [~,~,ext]=fileparts(filtermask);
    if strcmp(ext,'.gz')
        tmpfn=tempname;
        copyfile(filtermask,[tmpfn,'.nii.gz']);
        gunzip([tmpfn,'.nii.gz']);
        ea_delete([tmpfn,'.nii.gz']);
        filtermask=[tmpfn,'.nii'];
    end
    filter=ea_open_vol(filtermask); 
    fibvox=fibers;
    fibvox(:,4)=1;
    fibvox=(filter.mat\fibvox');
    ids=spm_sample_vol(filter,double(fibvox(1,:)),double(fibvox(2,:)),double(fibvox(3,:)),1);
    fibers=fibers(logical(ids),:);
    
    % recompute idx:
    [~,iax,iac]=unique(fibers(:,4));
    clear idx
    for fib=1:length(iax)-1
        idx(fib,1)=iax(fib+1)-iax(fib);
    end
    % add last entry
    idx(fib+1,1)=sum(fibers(:,4)==max(fibers(:,4)));
    fibers(:,4)=iac; % make sure no single fiber was completely erased and is not accounted for.
end

% Use other meta info from the first FTR file
ftr = load(ftrFiles{1});
ftr.fibers = fibers;
ftr.idx = idx;
save(outputFile, '-struct', 'ftr', '-v7.3');
ea_ftr2trk(outputFile);
