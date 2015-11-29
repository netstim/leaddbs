function bbstruct=ea_viz2brainbrowser(varargin)

% Export function to brainbrowser made at OHBM Hackathon 2014 together with Tarek
% Sherif.

vizstruct=varargin{1};
if nargin==2
   jtype=varargin{2};
   bbstruct.type=jtype;
end

bbstruct.vertices=[];
bbstruct.normals=[];

bbstruct.colors=[];
offset=-1; % since indices should start at 0.

vizstruct=vizstruct(:);

for entry=1:length(vizstruct)
    
    bbstruct.vertices=[bbstruct.vertices;vizstruct(entry).vertices];
    
    for n=1:length(vizstruct(entry).normals);
        vizstruct(entry).normals(n,:)=(vizstruct(entry).normals(n,:)/norm(vizstruct(entry).normals(n,:)));
    end
    
    bbstruct.normals=[bbstruct.normals;vizstruct(entry).normals];
    if ~isempty(vizstruct(entry).vertices)
        if length(vizstruct(entry).colors(1,:))==4
            bbstruct.colors=[bbstruct.colors;repmat([vizstruct(entry).colors(1,:)],length(vizstruct(entry).vertices),1)];
        elseif length(vizstruct(entry).colors(1,:))==3
            bbstruct.colors=[bbstruct.colors;repmat([vizstruct(entry).colors(1,:),0.7],length(vizstruct(entry).vertices),1)]; 
        end
    end
    
    bbstruct.shapes(entry).indices=vizstruct(entry).faces+offset;
    
    bbstruct.shapes(entry).indices=bbstruct.shapes(entry).indices;
    offset=offset+length(vizstruct(entry).vertices);
end

bbstruct.vertices=bbstruct.vertices';
bbstruct.vertices=bbstruct.vertices(:);
bbstruct.colors=bbstruct.colors';
bbstruct.colors=bbstruct.colors(:);
bbstruct.normals=bbstruct.normals';
bbstruct.normals=bbstruct.normals(:);
bbstruct.normals(isnan(bbstruct.normals))=0;