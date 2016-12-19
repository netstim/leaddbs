function ftrstruct_merge(filenames,selection,outname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ftrstruct_merge(cellfilenames, selections, ftroutname)
%
% merges selections of several ftrstructs into one ftrstruct-file
%
% params - cellfilenames contains a cellarray of ftrstruct-filenames, like
%          { 'folder1/patient1_FTR.mat', 'folder2/patient2_FTR.mat' }
%
%        - selections is an array of the same size as cellfilenames 
%          containing the fiber bundles as numbers, where 1 corresponds
%          to the whole tractgram, and numbers above one to selected
%          bundels, for example, selections = [1 2], constructs a ftrstruct
%          containing the whole tractogram of the first ftrstruct
%          and the first selected bundle of the second.
%
%        - ftroutname is the filename of the merged ftr
%
%
% example >> ftrstruct_merge( { 'folder1/patient1_FTR.mat', 'folder2/patient2_FTR.mat' },[1 2],'merged_FTR.mat');
%
% Marco Reisert 2012 UKLFR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





cnt = 1;
for k = 1:length(filenames),
    ftr = ftrstruct_read(filenames{k});  
    if k == 1,
        ftr_merged = ftr;
        ftr_merged.curveSegCell = [];
        ftr_merged.fiber = [];
    end;
    if selection(k) == 1,
        ftr_merged.curveSegCell = [ftr_merged.curveSegCell ftr.curveSegCell];
        nf = length(ftr.curveSegCell);
        ftr_merged.fiber{k}.name = filenames{k};
        ftr_merged.fiber{k}.curveID = cnt:(cnt+nf-1);
        ftr_merged.fiber{k}.roiName = {'fromfile'};
        ftr_merged.fiber{k}.user = [];        
        cnt = cnt + nf;
    else
        ftr_merged.curveSegCell = [ftr_merged.curveSegCell ftr.curveSegCell(ftr.fiber{selection(k)-1}.curveID)];
        nf = length(ftr.fiber{selection(k)-1}.curveID);
        ftr_merged.fiber{k}.name = filenames{k};
        ftr_merged.fiber{k}.curveID = cnt:(cnt+nf-1);
        ftr_merged.fiber{k}.roiName = {'fromfile'};
        ftr_merged.fiber{k}.user = [];        
        cnt = cnt + nf;                
    end;        
end;
 ftr_merged.connectCell = arrayfun(@(x) {x},1:length(ftr_merged.curveSegCell));
ftrstruct_write(ftr_merged,outname);
    