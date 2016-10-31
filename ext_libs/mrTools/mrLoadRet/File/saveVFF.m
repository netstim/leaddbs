% saveVFF.m
%
%        $Id$
%      usage: saveVFF(filename,data)
%         by: justin gardner
%       date: 10/19/07
%    purpose: save a VFF format (I don't actually know this format
%             backed it out of Jonas code so this is probably not general).
%
function retval = saveVFF(filename,data)

% check arguments
if ~any(nargin == [1 2])
  help saveVFF
  return
end

% add vff to filename
filename = sprintf('%s.vff',stripext(filename));

% squeeze the data
data = squeeze(data);
datalen = length(data);

% open the file
f = fopen(filename,'w','b');

% how to offset/scale values
maxValue=2^15-1;
scale = max(abs(data))/maxValue;
offset = maxValue;
data = data/scale+offset;

% write out header
fprintf(f,'ncaa\n');
fprintf(f,'rank=3;\n');
fprintf(f,'type=raster;\n');
fprintf(f,'format=slice;\n');
fprintf(f,'size=1 1 %i;\n',datalen);
fprintf(f,'origin=-1.000000 -1.000000 -0.000000;\n');
fprintf(f,'extent=0 0 %i;\n',datalen-1);
fprintf(f,'aspect=1.000000 1.000000 1.000000;\n');
fprintf(f,'bands=1;\n');
fprintf(f,'bits=16;\n');
fprintf(f,'value=%i %f;\n',offset,scale);
fprintf(f,'byte_order=big_endian;\n');
fprintf(f,'data_origin=left_anterior_superior;\n');
fprintf(f,'data_order=yzx;\n');
fwrite(f,12,'uint8');
fprintf(f,'\n');

% write out data
fwrite(f,data,'uint16');
fclose(f);

