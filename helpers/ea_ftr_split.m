function ea_ftr_split(ftrFile, numSubj, outputFolder, subjIDs)
% Split an aggregated connectome into subject-wise connectomes

arguments
    ftrFile         {mustBeFile}            % input FTR file
    numSubj         {mustBeNumeric}         % number of subjects
    outputFolder    {mustBeTextScalar} = '' % output folder
    subjIDs         {mustBeText} = {}       % ID of subjs
end

ftrFile = GetFullPath(ftrFile);

load(ftrFile, 'idx');
if mod(length(idx), numSubj)
    ea_error('Number of subjects is not correct!', showdlg=0, simpleStack=1);
end

if isempty(outputFolder)
    [ftrFolder, ftrName] = fileparts(ftrFile);
    outputFolder = fullfile(ftrFolder, [ftrName, '_SubjWise']);
    ea_mkdir(outputFolder);
end

fmtStr = ['%0', num2str(numel(num2str(numSubj))), 'd'];
if isempty(subjIDs)
    subjIDs = strcat('subj', num2str((1:numSubj)', fmtStr));
end

ea_cprintf('Comments', 'Loading connectome...\n');
ftr = load(ftrFile);

ea_cprintf('Comments', 'Splitting connectome...\n');

numFibPerSubj = length(ftr.idx)/numSubj;
startIndex = 1:numFibPerSubj:length(ftr.idx);
endIndex = startIndex + numFibPerSubj - 1;

subjFTR = struct('ea_fibformat', '1.1', 'fourindex', 1);
if isfield(ftr, 'voxmm')
    subjFTR.voxmm = ftr.voxmm;
end
if isfield(ftr, 'mat')
    subjFTR.mat = ftr.mat;
end

cumIdx = [0; cumsum(ftr.idx)];
for i=1:numSubj
    fibStartIndex = cumIdx(startIndex(i)) + 1;
    fibEndIndex = cumIdx(endIndex(i)+1);

    subjFTR.fibers = ftr.fibers(fibStartIndex:fibEndIndex,:);
    subjFTR.idx = ftr.idx(startIndex(i):endIndex(i));

    save(fullfile(outputFolder, subjIDs{i}), '-struct', 'subjFTR', '-v7.3');

    disp(['Subj ', num2str(i, fmtStr), '/', num2str(numSubj), '...']);
end
