function bbstruct=ea_viz2brainbrowser(vizstruct)

% Export function to brainbrowser made at OHBM Hackathon 2014 together with Tarek
% Sherif.

bbstruct.vertices=[];
bbstruct.normals=[];

bbstruct.colors=[];
offset=-1; % since indices should start at 0.

vizstruct=vizstruct(:);

for entry=1:length(vizstruct)
    
    bbstruct.vertices=[bbstruct.vertices;vizstruct(entry).vertices];
    
    for n=1:length(vizstruct(entry).normals);
        vizstruct(entry).normals(n,:)=(vizstruct(entry).normals(n,:)/norm(vizstruct(entry).normals(n,:)))*-1;
    end
    bbstruct.normals=[bbstruct.normals;vizstruct(entry).normals];
    if ~isempty(vizstruct(entry).vertices)
        bbstruct.colors=[bbstruct.colors;repmat([vizstruct(entry).colors(1,:)],length(vizstruct(entry).vertices),1)];
        
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