function cfv=ea_mapcolvert2face(varargin)

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

% check whether fvcdat is size of faces
if isequal(size(cfv.facevertexcdata),size(cfv.faces))
    return
else
    nufvcdat=zeros(size(cfv.faces));
    for face=1:size(cfv.faces,1)
        nufvcdat(face,:)=mean(cfv.facevertexcdata(cfv.faces(face,:),:));
    end
    cfv.facevertexcdata=nufvcdat;
    
end