function [data,hdr]=loadVFF( filename, loadOnlyHeader )
% [data,hdr]=tfiReadVFF( filename, <loadOnlyHeader>)
% 
% Loads data in .vff format
%
% this function is borrowed from Jonas' TFI distribution
% 

% default to load data and hdr
if ieNotDefined('loadOnlyHeader'),loadOnlyHeader=0;end


endofhdr=0;
hdr.byteorder='b';

data=[];hdr=[];
filename = sprintf('%s.vff',stripext(filename));
if ~isfile(filename)
  disp(sprintf('(loadVFF) Could not open %s',filename));
  return
end
f=fopen( filename, 'rb');
str=fgetl(f);

while (~endofhdr & ~isempty(str))
  [a,c]=sscanf(str,'size=%i %i %i;');
  if (c==3)
    hdr.size=[a(2) a(1) a(3)];
  end
  [a,c]=sscanf(str,'aspect=%f %f %f;');
  if (c==3)
    hdr.aspect=[a(2) a(1) a(3)];
  end
  [a,c]=sscanf(str,'value=%i %f;');
  if (c==2)
    hdr.scale=a(2);
    hdr.offset=a(1);
  end
  [a,c]=sscanf(str,'bits=%i;');
  if (c==1)
    hdr.bits=a(1);
    hdr.bytes=a(1)/8;
  end
  [a,c]=sscanf(str,'byte_order=%s');
  if (c==1)
    switch (a)
     case 'little_endian;'
      hdr.byteorder='l';
     case 'big_endian;'
      hdr.byteorder='b';
    end
  end
  if (length(str)==1 & double(str)==12)
    endofhdr=1;
    break
  end
  str=fgetl(f);
end
fclose(f);
if loadOnlyHeader,return,end

prec=['uint' num2str(hdr.bits)];

f=fopen( filename, 'rb',hdr.byteorder);
while (double(fgetl(f))~=12)
end

data=fread(f,prod(hdr.size),prec);
fclose(f);

% scale - only for 16-bit ushorts
if (hdr.bits==16)
  VFF_MAX_USHORT=2^16-1;
  data=(data-(VFF_MAX_USHORT-hdr.offset))*hdr.scale;
end

% rearrange into x,y,z order
data=permute(data,[2,1,3]);
