function [fName, path, errStr]= split_filename(longFileName)
% function [fName, path, errStr]= split_filename(longFileName)
%
%
%  Bjoern W. Kreher
%  03/04
% 
%  UNIX

fName= '';   path= '';   errStr= '';

if ispc
    idx= find(longFileName == '\');
elseif isunix
    idx= find(longFileName == '/');
else
    idx= find((longFileName == '/') | (longFileName == '\'));
end
if ~isempty(idx) && (idx(end)+1 < length(longFileName)) && isdir(longFileName(1:idx(end)))
    path= longFileName(1:idx(end));
    fName= longFileName((idx(end) + 1):end);
else
    fName= longFileName;
    errStr= 'Can''t determine directory tag';
end
