function [List,ListSize,ListDate] = LeG_findSubFiles(RootDir,SearchStr)

% Searching RootDir
% h = waitbar(0,'Searching...');
a = cell2mat(cellfun(@genpath,{RootDir},'uniformoutput',false));
if regexpi(computer,'win')
    b = regexp(a,';','split');
else
    b = regexp(a,':','split');
end
List = {};
ListSize = [];
ListDate = [];
for k=1:length(b)
%     if ~rem(k,10) %execute every 10 iterations
%         waitbar(k/length(b),h);
%     end
    c = dir(b{k});
    list = {c(~[c.isdir]).name};
    fsize = ([c(~[c.isdir]).bytes]./1e9)';
    fdate = [c(~[c.isdir]).datenum]';
    list = cellfun(@(x) fullfile(b{k},x),list,'uniformoutput',false)';
    List = [List;list];    
    ListSize = [ListSize;fsize];
    ListDate = [ListDate;fdate];
end
if ~isempty(SearchStr) %search a subset based on FileType
    SearchBoolean = ~cellfun(@isempty,regexp(List,SearchStr));
    List = List(SearchBoolean);
    ListSize = ListSize(SearchBoolean);
    ListDate = ListDate(SearchBoolean);
end
% close(h);
