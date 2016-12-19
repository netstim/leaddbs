% function erg = save_brooker(fName, mrStruct)
%
%  save_brooker supports no echos until now
%
%

function erg = save_brooker(fName, mrStruct)

% errors and tests
% if mrstruct_istype(mrStruct, 'series3D') ~= 1
%    warning('save_brooker supports no echos until now');
%    return;
% end


% preallocation of X
% end of preallocation of X

% read of the data

fid=fopen(fName,'wb','ieee-be');
fseek(fid,0,-1); 
dataAy= reshape(permute(mrStruct.dataAy, [2 1 3 4]), [prod(size(mrStruct.dataAy)) 1]);
fwrite(fid, dataAy, 'long');
clear dataAy
fclose(fid);
erg= fName;
