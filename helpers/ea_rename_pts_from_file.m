function ea_rename_pts_from_file(BIDSRoot, renameListFile)
% Function to rename subj in a BIDS dataset from a renaming list file

firstLine = fgetl(fopen(renameListFile, 'rt'));
if isempty(regexpi(firstLine, '^oldSubjId.*newSubjId$', 'once'))
    table = readtable(renameListFile, 'NumHeaderLines', 0);
else
    table = readtable(renameListFile, 'NumHeaderLines', 1);
end

subjID = table2cell(table);

ea_rename_pts(BIDSRoot, subjID(:,1), subjID(:,2));
