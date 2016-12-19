function [fileStr, errStr]= str2filename(nameStr, structFlag)
%function [fileStr, errStr]= str2filename(nameStr, structFlag)
%%   
%
% Bjoern W. Kreher
% 3/04
%
% UNIX

fileStr= '';errStr= '';

allowedChar= double('.0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ-_abcdefghijklmnopqrstuvwxyz');
if exist('structFlag') && ~isempty(structFlag)
    if strcmp(structFlag, 'extFileName')
        allowedChar= double('.0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ-_abcdefghijklmnopqrstuvwxyz ()!');
    else
        allowedChar= double('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz');
    end
end

[X, Y]=meshgrid(double(nameStr), double(allowedChar));
idx= sum(X == Y, 1) == 0;

fileStr= nameStr;
fileStr(idx)= '_';