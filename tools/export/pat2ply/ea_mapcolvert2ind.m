function cfv=ea_mapcolvert2ind(varargin)

cfv=varargin{1};
if nargin>1
    jetlist=varargin{2};
else % assume default color map
    jetlist=jet; 
end

% check if fvcdat is rgb
if size(cfv.facevertexcdata,2)==1
   cfv.facevertexcdata=squeeze(ind2rgb(round(cfv.facevertexcdata),jetlist));
end

% scale to uint8
cfv.facevertexcdata=uint8(cfv.facevertexcdata*255);