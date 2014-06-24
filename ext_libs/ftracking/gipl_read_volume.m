function [V] = gipl_read_volume(info)
% function for reading volume of Guys Image Processing Lab (Gipl) volume file
% 
% volume = gipl_read_header(file-header)
%
% examples:
% 1: info = gipl_read_header()
%    V = gipl_read_volume(info);
%    imshow(squeeze(V(:,:,round(end/2))),[]);
%
% 2: V = gipl_read_volume('test.gipl');

if(~isstruct(info)) info=gipl_read_header(info); end

% Open gipl file
f=fopen(getfield(info,'filename'),'rb','ieee-be');

  % Seek volume data start
  if(info.image_type==1), voxelbits=1; end
  if(info.image_type==7||info.image_type==8), voxelbits=8; end
  if(info.image_type==15||info.image_type==16), voxelbits=16; end
  if(info.image_type==31||info.image_type==32||info.image_type==64), voxelbits=32; end
  if(info.image_type==65), voxelbits=64; end
  
  datasize=prod(getfield(info,'sizes'))*(voxelbits/8);
  fsize=getfield(info,'filesize');
  fseek(f,fsize-datasize,'bof');

  % Read Volume data
  volsize(1:3)=getfield(info,'sizes');

  if(info.image_type==1), V = logical(fread(f,datasize,'bit1')); end
  if(info.image_type==7), V = int8(fread(f,datasize,'char')); end
  if(info.image_type==8), V = uint8(fread(f,datasize,'uchar')); end
  if(info.image_type==15), V = int16(fread(f,datasize,'short')); end 
  if(info.image_type==16), V = uint16(fread(f,datasize,'ushort')); end
  if(info.image_type==31), V = uint32(fread(f,datasize,'uint')); end
  if(info.image_type==32), V = int32(fread(f,datasize,'int')); end
  if(info.image_type==64), V = single(fread(f,datasize,'float')); end 
  if(info.image_type==65), V = double(fread(f,datasize,'double')); end 

fclose(f);

% Reshape the volume data to the right dimensions
V = reshape(V,volsize);

