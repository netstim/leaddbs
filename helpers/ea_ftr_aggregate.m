function [fibers, idx, partition] = ea_ftr_aggregate(ftrFiles, outputFile, opts)
% Aggregate fibers in the input FTR files
%
% Example:
%     ea_ftr_aggregate(ftrFiles, 'FTR.mat');                              % Aggregate fibers in ftrFiles
%     ea_ftr_aggregate(ftrFiles, 'FTR.mat', sel=20000, type='number');    % Sample 20000 fibers and aggregate
%     ea_ftr_aggregate(ftrFiles, 'FTR.mat', sel=2, type='interval');      % Sample with interval 2 and aggregate
%     ea_ftr_aggregate(ftrFiles, 'FTR.mat', sel=1:2:20000, type='index'); % Sample the specified fibers and aggregate

arguments
    ftrFiles        {mustBeText}                                                 % Input FTR files
    outputFile      {mustBeTextScalar}                                           % Output file name
    opts.sel        {mustBeInteger}                                              % number/interval/index of selected fibers
    opts.type       {mustBeMember(opts.type, {'number', 'interval', 'index'})}   % Type of 'opts.sel'
    opts.mask       {mustBeFile}                                                 % Filter mask for output
    opts.ref        {mustBeFile}                                                 % Reference space when exporting trk
end

ftrFiles = GetFullPath(ftrFiles);
if ~all(isfile(ftrFiles))
    ea_error(sprintf('Missing FTR file[s]:\n%s\n', strjoin(ftrFiles(~isfile(ftrFiles)), '\n')), showdlg=0, simpleStack=1);
end

fibers = [];
idx = [];
partition = []; % Store Num fibers of each FTR file

for i=1:length(ftrFiles)
    fprintf('Aggregating fibers from %d/%d FTR files...\n', i, length(ftrFiles));
    if ~isfield(opts, 'sel') % Simple aggregate
        ftr = load(ftrFiles{i});
    else % Sample and aggregate
        if isscalar(opts.sel) % sel is a single number
            if ~isfield(opts, 'type') % Treat as sample number by default if type not specified
                warning('Number type not specified! Used as sample number.')
                opts.type = 'number';
            end
        else % sel is a serial of numbers
            opts.type = 'index';
        end
        [ftr.fibers, ftr.idx] = ea_ftr_sample(ftrFiles{i}, opts.sel, opts.type);
    end

    fibers = [fibers;ftr.fibers];
    idx = [idx;ftr.idx]; 
    partition = [partition; length(ftr.idx)];
end

% Reset fibers(:, 4)
fibers(:, 4) = repelem(1:numel(idx), idx);

% filter
if isfield(opts, 'mask')
    [~,~,ext]=fileparts(opts.mask);
    if strcmp(ext,'.gz')
        tmpfn=tempname;
        copyfile(opts.mask,[tmpfn,'.nii.gz']);
        gunzip([tmpfn,'.nii.gz']);
        ea_delete([tmpfn,'.nii.gz']);
        opts.mask=[tmpfn,'.nii'];
    end
    filter=ea_open_vol(opts.mask); 
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

if ~isfield(opts, 'ref')
    spacedef = ea_getspacedef;
    opts.ref = [ea_space, spacedef.templates{1}, '.nii'];
end

ea_ftr2trk(outputFile, opts.ref);
